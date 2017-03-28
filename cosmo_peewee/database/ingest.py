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
from .models import Files, NUV_raw_headers, NUV_corr_headers, FUV_primary_headers, FUVA_raw_headers, FUVB_raw_headers
from .models import FUVA_corr_headers, FUVB_corr_headers, Lampflash, Rawacqs, Darks, Stims, Observations

from .database_keys import nuv_raw_keys, nuv_corr_keys, fuv_primary_keys, fuva_raw_keys, fuvb_raw_keys
from .database_keys import fuva_corr_keys, fuvb_corr_keys, obs_keys, file_keys

from ..filesystem import find_all_datasets

from ..osm.monitor import pull_flashes
from ..osm.monitor import monitor as osm_monitor

from ..dark.monitor import monitor as dark_monitor
from ..dark.monitor import pull_orbital_info

#-------------------------------------------------------------------------------

def bulk_insert(table, data_source):
    
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


    #-- Convert itertools.chain object to list...
    data_source = list(data_source)

    database = get_database()
    database.connect()

    try:
        with database.atomic():
            #-- Only add 100 files at a time....
            for idx in range(0, len(data_source), 100):
                table.insert_many(data_source[idx:idx+100]).execute()

    #-- Lots of multiples, hopefully will be fixed with new filesystem implimentation.   
    except IntegrityError as e:
        print('IntegrityError:', e)
    except IOError as e:
        print('IOError:',  e)
    except OperationalError as e:
        print('OperationalError', e)
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
        
        #-- Some files exist in files that need to be removed that don't exist in Observations.
        try:
            Observations.get(Observations.filename == f).delete_instance()
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
        
    if len(data_to_insert):
        #-- Pass to bulk insert.    
        bulk_insert(Files, itertools.chain(*data_to_insert))

#-------------------------------------------------------------------------------

def populate_observations(num_cpu=2):
    
    """Populate table with observations. This seperates files that the 
    monitors use from telemetry/MAST created files like .jit, .jif, cci's, etc.

    Parameters
    ----------
    num_cpu: int
        number of worker processes.

    Returns
    -------
    None
    
    """

    search_list = ['%rawtag%.gz',
                   '%corrtag%.gz',
                   '%_x1d.%.gz', #-- Get x1ds
                   '%_x1ds%.gz', #-- Get x1dsums
                   '%rawacq%.gz%',
                   '%lampflash%.gz']
        
    for search_str in search_list:
        database = get_database()
        database.connect()

        logger.info('POPULATING {} INTO OBSERVATIONS'.format(search_str))
        files_to_add = (Files.select()
                            .where(
                                Files.filename.contains(search_str) & 
                                Files.filename.not_in(Observations.select(Observations.filename)) &
                                (Files.monitor_flag == True)
                            ))

        database.close() 
        
        #-- set up partials
        partial = functools.partial(pull_data,
                                    function=obs_keys)
        
        #-- pool up the partials and pass it the iterable (files)
        pool = mp.Pool(processes=num_cpu)
        data_to_insert = pool.map(partial, files_to_add)
        
        if len(data_to_insert):
            #-- Pass to bulk insert.
            bulk_insert(Observations, itertools.chain(*data_to_insert))
        
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
    
    #-- Have to do some special parsing for FUV because of the degeneracy in A & B 
    #-- data products. Dont want to include duplicates.
    if search_str in ['%rawtag_a.fits.gz%', '%rawtag_b.fits.gz%'] and table._meta.db_table == 'fuv_primary_headers':
        files_to_add = (Observations.select()
                        .where(
                            Observations.filename.contains(search_str) & 
                            Observations.rootname.not_in(table.select(table.rootname))
                        ))

    #-- Because the FUV and NUV x1d files share the same naming scheme, we need to 
    #-- perform a join on the resepective table and make sure that we are only adding
    #-- FUV or NUV exposures to their respective tables.
    elif table._meta.db_table == 'fuv_x1d_headers':
        files_to_add = (Observations.select()
                         .join(FUV_raw_headers)
                         .where(
                             Observations.filename.contains(search_str) & 
                             Observations.rootname.not_in(table.select(table.rootname)) &
                             (FUV_raw_headers.detector == 'FUV')
                         ))

    #-- NUV x1d's... Maybe try to combine with the FUV by using a dictionary.
    elif table._meta.db_table == 'nuv_x1d_headers':
        files_to_add = (Observations.select()
                         .join(NUV_raw_headers)
                         .where(
                             Observations.filename.contains(search_str) & 
                             Observations.rootname.not_in(table.select(table.rootname)) &
                             (NUV_raw_headers.detector == 'NUV') 
                         ))
    
    #-- Else, the tables should look like this
    else:
        files_to_add = (Observations.select()
                            .where(
                                Observations.filename.contains(search_str) & 
                                Observations.filename.not_in(table.select(table.filename)) 
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

    files_to_add = (Observations.select()
                            .where(
                                Observations.filename.contains('%lampflash%.gz') & 
                                Observations.filename.not_in(Lampflash.select(Lampflash.filename)) 
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

    files_to_add = (Observations.select()
                            .where(
                                Observations.filename.contains('%rawacq%.gz') & 
                                Observations.filename.not_in(Rawacqs.select(Rawacqs.filename))
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

    files_to_add = (Observations.select()
                            .where( 
                                Observations.filename.contains('%corrtag%.gz') & 
                                Observations.filename.not_in(Darks.select(Darks.filename)) & 
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

    files_to_add = (Observations.select()
                            .where(
                                Observations.filename.contains('%corrtag\_%.gz') & 
                                Observations.filename.not_in(Stims.select(Stims.filename)) 
                            ))
    database.close()

    partial = functools.partial(bulk_insert,
                                function=locate_stims)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    data_to_insert = pool.map(partial, files_to_add)
    
    if len(data_to_insert): 
        #-- Pass to bulk insert.
        bulk_insert(Stims, itertools.chain(*data_to_insert))

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
              NUV_raw_headers,
              NUV_corr_headers,
              FUV_primary_headers, 
              FUVA_raw_headers, 
              FUVB_raw_headers,
              FUVA_corr_headers, 
              FUVB_corr_headers,
              Lampflash,
              Rawacqs,
              Darks]
        
    #-- Safe checks existance of tables first to make sure they dont get clobbered.
    database.create_tables(tables, safe=True)
    
    #-- Close DB connection.
    database.close()
    
    #-- Files
    logger.info("Ingesting Files from {}".format(settings['data_location']))
    populate_files(settings)

    #-- Observation table    
    # logger.info("Populating observations table.")
    # populate_observations(settings['num_cpu'])

    # #-- NUV rawtag headers    
    # logger.info("Populating NUV rawtag headers.")
    # populate_tables(NUV_raw_headers, nuv_raw_keys, '%_rawtag.fits.gz%', settings['num_cpu'])

    # #-- NUV corrtag headers    
    # logger.info("Populating NUV corrtag headers.")
    # opulate_tables(NUV_corr_headers, nuv_corr_keys, '%_corrtag.fits.gz%', settings['num_cpu'])

    # #-- FUV primary headers    
    # logger.info("Populating FUV primary headers.")
    # populate_tables(FUV_primary_headers, fuv_primary_keys, '%rawtag_a.fits.gz%', settings['num_cpu'])
    # populate_tables(FUV_primary_headers, fuv_primary_keys, '%rawtag_b.fits.gz%', settings['num_cpu'])

    # #-- FUV rawtag headers    
    # logger.info("Populating FUV rawtag headers.")
    # populate_tables(FUVA_raw_headers, fuva_raw_keys, '%rawtag_a.fits.gz%', settings['num_cpu'])
    # populate_tables(FUVB_raw_headers, fuvb_raw_keys, '%rawtag_b.fits.gz%', settings['num_cpu'])

    # #-- FUV corrtag headers    
    # logger.info("Populating FUV corrtag headers.")
    # populate_tables(FUVA_corr_headers, fuva_corr_keys, '%corrtag_a.fits.gz%', settings['num_cpu'])
    # populate_tables(FUVB_corr_headers, fuvb_corr_keys, '%corrtag_b.fits.gz%', settings['num_cpu'])

    # #-- Populate rawacq monitor meta
    # logger.info("Populating Rawacqs Table")
    # populate_acqs(settings['num_cpu'])

    # #-- Populate OSM monitor metadata
    # logger.info("Populating OSM Shift Table")
    # populate_osm(settings['num_cpu'])

    # #-- Populate Darks monitor meta
    # logger.info("Populating Darks Table")
    # populate_darks(settings['num_cpu'])
    
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
    
    osm_monitor()
    dark_monitor()

#-------------------------------------------------------------------------------
