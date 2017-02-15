from __future__ import print_function, absolute_import, division

import os
import sys

import logging
logger = logging.getLogger(__name__)

from astropy.io import fits
from peewee import *
import functools
import multiprocessing as mp

from .models import get_database, get_settings
from .models import Files, NUV_raw_headers, NUV_corr_headers, FUV_primary_headers, FUVA_raw_headers, FUVB_raw_headers, FUVA_corr_headers, FUVB_corr_headers
from .database_keys import nuv_raw_keys, nuv_corr_keys, fuv_primary_keys, fuva_raw_keys, fuvb_raw_keys, fuva_corr_keys, fuvb_corr_keys
from ..filesystem import find_all_datasets

from ..osm.monitor import pull_flashes
from ..osm.monitor import monitor as osm_monitor

#-------------------------------------------------------------------------------

def bulk_insert(file_result, table, function):
    
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
        table.create(**data)
    except IntegrityError as e:
        print(e)
        print(os.path.join(file_result.path, file_result.filename))
    except IOError as e:
        print(e)
        print(os.path.join(file_result.path, file_result.filename))

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

    #-- Not the prettiest thing in the whole world....
    files_to_remove = [os.path.basename(full_path) for full_path in all_files if not os.path.exists(full_path)]
    
    database = get_database()
    database.connect()

    Files.delete().where(Files.filename.in_(files_to_remove)).execute()  
    
    database.close()      

#-------------------------------------------------------------------------------

def setup_logging():
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

    database = get_database()
    database.connect()
    
    #-- Get previous files
    previous_files = {os.path.join(result.path, result.filename) for result in Files.select(Files.path, Files.filename)}
    
    #-- Make sure that the previous files in the database actually exist.
    bulk_delete(previous_files)
    
    for (path, filename) in find_all_datasets(settings['data_location'], settings['num_cpu']):
        #-- Join path + files together
        full_filepath = os.path.join(path, filename)
       
        #-- If file is in table then move to the next iteration. If not, add.
        if full_filepath in previous_files:
            continue

        rootname = filename.split('_')[0]
        #-- Rootname for files are 9 characters long, if not the NULL
        if not len(rootname) == 9:
            rootname = None
        
        #-- Set monitoring flags for ingestion, we have some corrupt files in the filesystem
        #-- False flags will allow us to avoid having NULL rows in monitoring metadata tables
        #-- and also track them.
        monitor_flag = True
        try:
            fits.open(full_filepath)
        except IOError:
            monitor_flag = False

        #-- There are some fishy files in the filesystem we need to fish out, here is way to expose them.
        try:
            #-- Insert metadata
            Files.insert(path=path,
                         filename=filename,
                         rootname=rootname,
                         monitor_flag=monitor_flag).execute() 
        except IntegrityError:
            print(full_filepath)
    
    database.close()

#-------------------------------------------------------------------------------

def populate_tables(table, table_keys, search_str, num_cpu=2):

    #-- Setting up logging
    setup_logging()
    
    #-- Connect to DB
    database = get_database()
    database.connect()
    
    #-- Select files where the filename looks like %search_string% and filename is
    #-- not in the table you wish to populate and that the monitoring flag is set 
    #-- to true so we know that it isnt corrupted.
    
    #-- Have to do some special parsing for FUV because of the degeneracy in A & B 
    #-- data products. Dont want to include duplicates.
    if search_str in ['%rawtag_a.fits.gz%', '%rawtag_b.fits.gz%'] and table._meta.db_table == 'fuv_primary_headers':
        files_to_add = (Files.select()
                        .where(
                            Files.filename.contains(search_str) & 
                            Files.rootname.not_in(table.select(table.rootname)) &
                            (Files.monitor_flag == True)
                        ))

    #-- Because the FUV and NUV x1d files share the same naming scheme, we need to 
    #-- perform a join on the resepective table and make sure that we are only adding
    #-- FUV or NUV exposures to their respective tables.
    elif table._meta.db_table == 'fuv_x1d_headers':
        files_to_add = (Files.select()
                         .join(FUV_raw_headers)
                         .where(
                             Files.filename.contains(search_str) & 
                             Files.rootname.not_in(table.select(table.rootname)) &
                             (FUV_raw_headers.detector == 'FUV') &
                             (Files.monitor_flag == True)
                         ))

    #-- NUV x1d's... Maybe try to combine with the FUV by using a dictionary.
    elif table._meta.db_table == 'nuv_x1d_headers':
        files_to_add = (Files.select()
                         .join(NUV_raw_headers)
                         .where(
                             Files.filename.contains(search_str) & 
                             Files.rootname.not_in(table.select(table.rootname)) &
                             (NUV_raw_headers.detector == 'NUV') &
                             (Files.monitor_flag == True)
                         ))
    
    #-- Else, the tables should look like this
    else:
        files_to_add = (Files.select()
                            .where(
                                Files.filename.contains(search_str) & 
                                Files.filename.not_in(table.select(table.filename)) &
                                (Files.monitor_flag == True)
                            ))
    
    #-- Show the user how many files are being added to the DB tables.
    #logger.info('ADDING {} FILES TO {}'.format(files_to_add.count(), table._meta.db_table))
    
    #-- set up partials
    partial = functools.partial(bulk_insert,
                                table=table,
                                function=table_keys)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    pool.map(partial, files_to_add)
    
    database.close()

#-------------------------------------------------------------------------------
def populate_osm(**kwargs):
    
    database = get_database()
    database.connect()

    files_to_add = (Files.select()
                            .where(
                                (Files.filename.contains('%lampflash%') | Files.filename.contains('%_rawacq%')) & 
                                Files.filename.not_in(table.select(table.filename)) &
                                (Files.monitor_flag == True)
                            ))

    '''
    partial = functools.partial(bulk_insert,
                                table=table,
                                function=function)
    
    #-- pool up the partials and pass it the iterable (files)
    pool = mp.Pool(processes=num_cpu)
    pool.map(partial, files_to_add)
    '''
    database.close()
#-------------------------------------------------------------------------------

def ingest_all():
    """Create tables and run all ingestion scripts"""
    setup_logging()
    #-- Get settings for functions
    settings = get_settings()
    
    #-- Open DB
    database = get_database()
    database.connect()
    
    #-- Tables to create.
    tables = [Files,
              NUV_raw_headers,
              NUV_corr_headers,
              FUV_primary_headers, 
              FUVA_raw_headers, 
              FUVB_raw_headers,
              FUVA_corr_headers, 
              FUVB_corr_headers]


    #-- Safe checks existance of tables first to make sure they dont get clobbered.
    database.create_tables(tables, safe=True)
    
    #-- Close DB
    database.close()
    
    #-- Files
    logger.info("Ingesting Files from {}".format(settings['data_location']))
    populate_files(settings)

    #-- NUV rawtag headers    
    logger.info("Populating NUV rawtag headers.")
    populate_tables(NUV_raw_headers, nuv_raw_keys, '%_rawtag.fits.gz%', settings['num_cpu'])

    #-- NUV corrtag headers    
    logger.info("Populating NUV corrtag headers.")
    populate_tables(NUV_corr_headers, nuv_corr_keys, '%_corrtag.fits.gz%', settings['num_cpu'])

    #-- FUV primary headers    
    logger.info("Populating FUV primary headers.")
    populate_tables(FUV_primary_headers, fuv_primary_keys, '%rawtag_a.fits.gz%', settings['num_cpu'])
    populate_tables(FUV_primary_headers, fuv_primary_keys, '%rawtag_b.fits.gz%', settings['num_cpu'])

    #-- FUV rawtag headers    
    logger.info("Populating FUV rawtag headers.")
    populate_tables(FUVA_raw_headers, fuva_raw_keys, '%rawtag_a.fits.gz%', settings['num_cpu'])
    populate_tables(FUVB_raw_headers, fuvb_raw_keys, '%rawtag_b.fits.gz%', settings['num_cpu'])

    #-- FUV corrtag headers    
    logger.info("Populating FUV corrtag headers.")
    populate_tables(FUVA_corr_headers, fuva_corr_keys, '%corrtag_a.fits.gz%', settings['num_cpu'])
    populate_tables(FUVB_corr_headers, fuvb_corr_keys, '%corrtag_b.fits.gz%', settings['num_cpu'])

    #-- Populate OSM monitor meta
    logger.info("Populating OSM Shift Table")
    populate_osm(OSM_shift, pull_flashes, settings['num_cpu'])
#-------------------------------------------------------------------------------

