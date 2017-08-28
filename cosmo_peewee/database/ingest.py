from __future__ import print_function, absolute_import, division

import os
import sys
import types

import logging
logger = logging.getLogger(__name__)

from astropy.io import fits

import argparse

from peewee import *

import functools

import itertools

import multiprocessing as mp
import numpy as np
from .models import get_database, get_settings
from .models import Files, NUV_corr_headers, FUVA_raw_headers, FUVB_raw_headers
from .models import FUVA_corr_headers, FUVB_corr_headers, Lampflash, Rawacqs, Darks, Stims, Observations, Gain

from .database_keys import nuv_corr_keys, fuva_raw_keys, fuvb_raw_keys
from .database_keys import fuva_corr_keys, fuvb_corr_keys, obs_keys, file_keys

from ..filesystem import find_all_datasets

from ..osm.monitor import pull_flashes
from ..osm.monitor import monitor as osm_monitor

from ..dark.monitor import monitor as dark_monitor
from ..dark.monitor import pull_orbital_info

from ..stims.monitor import locate_stims
from ..stims.monitor import stim_monitor

from ..cci.gainmap import write_and_pull_gainmap
#-------------------------------------------------------------------------------

def bulk_insert(table, data_source, debug=False):
    
    """ Ingest data into database

    Parameters
    ----------
    table: Peewee table object
        Table that is going to be populated.
    data_source: list
        A list full of dictionaries to be ingested.

    Returns
    -------
    None
    
    """


    database = get_database()
    database.connect()

    #-- Insert all of the entries one by one for debugging purposes.
    if debug:
        for item in data_source:
            try:
                print(item)
                table.insert(**item).execute()
            except IntegrityError as e:
                print('IntegrityError:', e)
                print(item['filename'])
        database.close()
    
    #-- Bulk inserts.
    else:
        try:
            with database.atomic():
                table.insert_many(data_source).execute()
            database.close()
        #-- Lots of multiples, hopefully will be fixed with new filesystem implimentation.   
        except IntegrityError as e:
            print('IntegrityError:', e)
            database.close()
        except IOError as e:
            print('IOError:',  e)
            database.close()
        except OperationalError as e:
            print('OperationalError', e)
            database.close()
        except InternalError as e:
            print('InternalError', e)
            database.close()
#-------------------------------------------------------------------------------

def pull_data(file_result, function):
    
    """This function inserts the rows into the data base. We use the partials
    to add data quick using pool.map().

    Parameters
    ----------
    table : peewee table object
        DB table you wish to ingest data into.
    function : Function
        Function that extracts the metadata to insert
    files : List
        List of files you wish to extract metadata from.
    
    """

    try:
    
        data = function(file_result)

        if isinstance(data, dict):
            data = [data]
        elif isinstance(data, types.GeneratorType):
            data = list(data)
        
        return data
    
    except IOError as e:
        print('IOError:',  e)
    
#-------------------------------------------------------------------------------

def bulk_delete(all_files):
    
    """
    This function is a temporary fix for /smov/cos/Data/. Since the 'date-string'
    in the file path changes, the database will try to ingest the same files. The
    paths changes due to newly downlinked data or when data is reprocessed using
    new reference files.

    Parameters
    ----------
    all_files : list
        A list of the current files in the database.

    Returns
    -------
    None
   
    """

    #-- Combine tuples of path and filenames to check for existance....
    combined_paths = [os.path.join(path, filename) for (path, filename) in all_files]

    #-- Not the prettiest thing in the whole world....
    files_to_remove = [os.path.basename(full_path) for full_path in combined_paths if not os.path.exists(full_path)]
    
    database = get_database()
    database.connect()

    logger.info('FOUND {} FILES TO DELETE!'.format(len(files_to_remove)))
    #-- Crude implimentation but right now it works.
    #-- Delete instances of files in DB that don't exist in files and observations tables.
    #-- This makes sure that the most up to date files in DB are ingested.    
    for f in files_to_remove:
        logger.info("DELETING INSTANCES OF {}".format(f))
        
        Files.get(Files.filename == f).delete_instance()
        
        #-- Hack to use the rootname to delete.
        try:
            Observations.get(Observations.rootname == f[:9]).delete_instance()
        except DoesNotExist:
            continue

    database.close()      

#-------------------------------------------------------------------------------

def setup_logging():
    
    """
    Set up logging....

    Parameters
    ----------
    None

    Returns
    -------
    None

    """
    
    #-- create the logging file handler
    logging.basicConfig(filename="cosmos_monitors.log",
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        level=logging.DEBUG)

    #-- handler for STDOUT
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logging.getLogger().addHandler(ch)

#-------------------------------------------------------------------------------

def populate_files(settings):
    
    """
    Populate paths and files from filesystem.
    
    Parameters
    ----------
    settings : dictionary
        A dictionary of credentials for database and monitors
    
    Returns
    -------
    None
    
    """

    #-- Open DB and get previous files.
    database = get_database()
    database.connect()

    previous_files = {(result.path, result.filename) for result in Files.select(Files.path, Files.filename)}

    database.close()

    #-- Bulk delete files that arent in /smov/cos/Data....
    bulk_delete(previous_files)

    #-- Create a set of all the data sets in /smov/cos/Data/
    smov_cos_data = set(find_all_datasets(settings['data_location'], settings['num_cpu']))

    #-- Turn the difference of the sets into a list to pass to bulk_insert
    files_to_add = list(smov_cos_data - previous_files)
    
    logger.info('FOUND {} FILES TO ADD!'.format(len(files_to_add)))

    #-- set up partials
    partial = functools.partial(pull_data,
                                function=file_keys)
        
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=settings['num_cpu'])
    data_to_insert = pool.map(partial, files_to_add)    
 
    #-- DB times out if all tables are nuked.
    #-- Seems to be the best method but don't know how to get it to work...
    if len(data_to_insert):
        #-- Little if else logic to avoid Integrity Errors not allowing good files to be ingested
        if len(files_to_add) < 6000:
            step = 1
        else:
            step = 100

        #-- Pass to bulk insert.
        for idx in range(0, len(list(data_to_insert)), step):    
            bulk_insert(Files, itertools.chain(*data_to_insert[idx:idx+step]))
        
#-------------------------------------------------------------------------------

def populate_tables(table, table_keys, search_str, num_cpu=2):
    
    """Parse data from different modes to proper tables. 

    Parameters
    ----------
    table: Peewee table object
        Table you wish to populate
    table_keys: function
        Keywords that will populate fields in tables.
    search_str: str
        A SQL formatted search string
    num_cpu: int
        number of worker processes.

    Returns
    -------
    None

    """

    #-- Connect to DB
    database = get_database()
    database.connect()
    
    #-- Select files where the filename looks like %search_string% and filename is
    #-- not in the table you wish to populate and that the monitoring flag is set 
    #-- to true so we know that it isnt corrupted.
    
    #-- Because the FUV and NUV x1d files share the same naming scheme, we need to 
    #-- perform a join on the resepective table and make sure that we are only adding
    #-- FUV or NUV exposures to their respective tables.
    if table._meta.db_table == 'fuv_x1d_headers':
        files_to_add = (Observations.select()
                         .where(
                             Observations.filename.contains(search_str) & 
                             Observations.rootname.not_in(table.select(table.rootname)) &
                             (Observations.detector == 'FUV')
                         ))

    #-- NUV x1d's... Maybe try to combine with the FUV by using a dictionary.
    elif table._meta.db_table == 'nuv_x1d_headers':
        files_to_add = (Files.select()
                         .where(
                             Files.filename.contains(search_str) & 
                             Files.rootname.not_in(table.select(table.rootname)) &
                             (Observations.detector == 'NUV') 
                         ))
    #-- Else, the tables should look like this
    else:
        files_to_add = (Files.select()
                            .where(
                                Files.filename.contains(search_str) & 
                                Files.rootname.not_in(table.select(table.rootname)) &
                                (Files.monitor_flag == True)
                            ))
    database.close()
    
    #-- set up partials
    partial = functools.partial(pull_data,
                                function=table_keys)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    data_to_insert = pool.map(partial, files_to_add)
    
    if len(data_to_insert):
        #-- Pass to bulk insert.
        bulk_insert(table, itertools.chain(*data_to_insert))

#-------------------------------------------------------------------------------
def populate_osm(num_cpu=2):
    
    """ Populate the OSM table
    
    Parameters
    ----------
    num_cpu: int
        number of worker processes.

    Returns
    -------
    None

    """
    
    database = get_database()
    database.connect()

    files_to_add = (Files.select()
                            .where(
                                Files.filename.contains('%lampflash%.gz') & 
                                Files.filename.not_in(Lampflash.select(Lampflash.filename)) & 
                                (Files.monitor_flag == True) 
                            ))
    database.close()

    partial = functools.partial(pull_data,
                                function=pull_flashes)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    data_to_insert = pool.map(partial, files_to_add)
    
    if len(data_to_insert):
        #-- Pass to bulk insert.
        bulk_insert(Lampflash, itertools.chain(*data_to_insert))

#-------------------------------------------------------------------------------
def populate_acqs(num_cpu=2):
    
    """ Populate the rawacqs table
    
    Parameters
    ----------
    num_cpu: int
        number of worker processes.

    Returns
    -------
    None

    """
    
    database = get_database()
    database.connect()

    files_to_add = (Files.select()
                            .where(
                                Files.filename.contains('%rawacq%.gz') & 
                                Files.filename.not_in(Rawacqs.select(Rawacqs.filename)) &
                                (Files.monitor_flag == True)
                            ))
    
    database.close()
 
    partial = functools.partial(pull_data,
                                function=pull_flashes)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    data_to_insert = pool.map(partial, files_to_add)
    
    if len(data_to_insert):
        #-- Pass to bulk insert.
        bulk_insert(Rawacqs, itertools.chain(*data_to_insert))

#-------------------------------------------------------------------------------

def populate_darks(num_cpu=2):
    
    """ Populate the dark table
    
    Parameters
    ----------
    num_cpu: int
        number of worker processes.

    Returns
    -------
    None

    """
    
    database = get_database()
    database.connect()

    files_to_add = (Files.select()
                            .join(Observations, on=(Files.rootname == Observations.rootname))
                            .where(
                                (Files.filename.contains('%corrtag%.gz') &
                                 Files.filename.not_in(Darks.select(Darks.filename)) & 
                                 Files.monitor_flag == True) &
                                (Observations.targname == 'DARK')
                            ))
    database.close()

    partial = functools.partial(pull_data,
                                function=pull_orbital_info)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    data_to_insert = pool.map(partial, files_to_add)

    if len(data_to_insert): 
        #-- Pass to bulk insert.
        bulk_insert(Darks, itertools.chain(*data_to_insert))

#-------------------------------------------------------------------------------

def populate_stims(num_cpu=2):
    
    """ Populate the stim pulse table
    
    None of the monitoring code has been integrated into the system yet.
    Still working with the monitoring team to refactor monitor....
    
    Parameters
    ----------
    num_cpu: int
        number of worker processes.

    Returns
    -------
    None

    """
    
    database = get_database()
    database.connect()

    files_to_add = (Files.select()
                            .where(
                                Files.filename.contains('%corrtag\_%.gz') & 
                                Files.filename.not_in(Stims.select(Stims.filename)) &
                                Files.monitor_flag == True
                            ))
    database.close()
    
    partial = functools.partial(pull_data,
                                function=locate_stims)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    
    #-- Add 100 files at a time...
    step=100
    for idx in range(0, len(list(files_to_add)), step):    
        data_to_insert = pool.map(partial, files_to_add[idx:idx+step])
        
        if len(data_to_insert):
            #-- Pass to bulk insert.
            bulk_insert(Stims, itertools.chain(*data_to_insert))
#-------------------------------------------------------------------------------

def populate_gain(num_cpu=2):
    
    """ 
    Populate gain table and create gainmaps
    
    Parameters
    ----------
    num_cpu: int
        number of worker processes.

    Returns
    -------
    None

    """
    
    database = get_database()
    database.connect()

    files_to_add = (Files.select()
                            .where(
                                (Files.filename.contains('l\_%\_00\____\_cci.fits%') |
                                Files.filename.contains('l\_%\_01\____\_cci.fits%')) & 
                                Files.filename.not_in(Gain.select(Gain.filename)) &
                                Files.monitor_flag == True
                            ))

    database.close()
    
    partial = functools.partial(pull_data,
                                function=write_and_pull_gainmap)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    print(len(files_to_add))
    #-- Add 10 files at a time.
    step=2
    for idx in range(0, len(list(files_to_add)), step):
        data_to_insert = pool.map(partial, files_to_add[idx:idx+step])
        bulk_insert(Gain, itertools.chain(*data_to_insert), debug=True)
        
    #     if len(data_to_insert):
    #         #-- Pass to bulk insert.
    #         bulk_insert(Gain, itertools.chain(*data_to_insert), debug=True)

    #-- For testing
    # files_to_add = Files.select().where(Files.filename.contains('l_2009226105625_01_167_cci.fits.gz'))
    # database.close()
    # data_to_insert = pool.map(partial, files_to_add)
    # bulk_insert(Gain, itertools.chain(*data_to_insert))
    #-- End testing
#-------------------------------------------------------------------------------

def ingest_all():
   
    """Create tables and run all ingestion scripts
    
    Parameters
    ----------
    None

    Returns
    -------
    None

    """
    
    setup_logging()
    
    #-- Get settings for functions
    settings = get_settings()
    
    #-- Open DB connection.
    database = get_database()
    database.connect()
    
    #-- Tables to create.
    tables = [Files,
              Observations,
              NUV_corr_headers,
              FUVA_raw_headers, 
              FUVB_raw_headers,
              FUVA_corr_headers, 
              FUVB_corr_headers,
              Lampflash,
              Rawacqs,
              Darks,
              Stims,
              Gain]

    #-- Safe checks existance of tables first to make sure they dont get clobbered.
    database.create_tables(tables, safe=True)
    
    #-- Close DB connection.
    database.close()
    
    #-- Files
    logger.info("INGESTING FILES FROM: {}".format(settings['data_location']))
    populate_files(settings)

    #-- Observation table    
    #-- Replacing the fuv_primary_headers and nuv_raw_headers and combining.
    logger.info("POPULATING OBSERVATION TABLE")

    #-- if raw file types don't exist, then whats the point...?
    filetypes = ['%rawacq.fits.gz',
                 '%lampflash.fits.gz',
                 '%rawaccum.fits.gz',
                 '%rawtag.fits.gz', 
                 '%rawtag_a.fits.gz', 
                 '%rawtag_b.fits.gz'
                 ]

    for file_type in filetypes:
        populate_tables(Observations, obs_keys, file_type, settings['num_cpu'])

    #-- NUV corrtag headers    
    logger.info("POPULATING NUV CORRTAGS")
    populate_tables(NUV_corr_headers, nuv_corr_keys, '%_corrtag.fits.gz%', settings['num_cpu'])

    #-- FUV rawtag headers    
    logger.info("POPULATING FUV RAWTAG HEADERS")
    populate_tables(FUVA_raw_headers, fuva_raw_keys, '%rawtag_a.fits.gz%', settings['num_cpu'])
    populate_tables(FUVB_raw_headers, fuvb_raw_keys, '%rawtag_b.fits.gz%', settings['num_cpu'])

    #-- FUV corrtag headers    
    logger.info("POPULATING FUV CORRTAGS")
    populate_tables(FUVA_corr_headers, fuva_corr_keys, '%corrtag_a.fits.gz%', settings['num_cpu'])
    populate_tables(FUVB_corr_headers, fuvb_corr_keys, '%corrtag_b.fits.gz%', settings['num_cpu'])

    #-- Populate rawacq monitor meta
    logger.info("POPULATING RAW ACQS TABLE")
    populate_acqs(settings['num_cpu'])

    #-- Populate OSM monitor metadata
    logger.info("POPULATING OSM DRIFT TABLE")
    populate_osm(settings['num_cpu'])

    #-- Populate Darks monitor meta
    logger.info("POPULATING DARKS TABLE")
    populate_darks(settings['num_cpu'])
    
    # #-- Populate Stim monitor
    # logger.info("POPULATING STIM TABLE")
    # #populate_stims(settings['num_cpu'])

    #-- Populate gain monitor
    logger.info("POPULATING GAIN TABLE")
    populate_gain(settings['num_cpu'])

#-------------------------------------------------------------------------------

def run_monitors():
    """ Run all COS Monitors
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    
    """
    setup_logging()
    
    #osm_monitor()
    dark_monitor()
    #stim_monitor()

#-------------------------------------------------------------------------------
