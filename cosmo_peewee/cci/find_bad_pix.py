from __future__ import absolute_import, print_function

"""
"""

__author__ = 'Mees Fix, Justin Ely'
__maintainer__ = 'Mees Fix'
__email__ = 'mfix@stsci.edu'
__status__ = 'Active'

import os
import glob
import logging
logger = logging.getLogger(__name__)

from peewee import *

from astropy.io import fits
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import multiprocessing as mp

import scipy
from scipy.optimize import leastsq, newton, curve_fit

from .constants import Y_BINNING, X_BINNING

from ..database.models import get_database, get_settings, Files
from ..database.models import Flagged_Pixels, Gain

import collections
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
def find_flagged(args):
    
    """
    Search the Gain table in cosmo for sagged pixels. Once located, there can be
    multiple entries because the gain for a pixel can drop below 3. We want the 
    date that it went bad first, this will become an entry in the Flagged_Pixel table.
    
    Parameters
    ----------
    args: Dict
        Dictionary containing HV_LVL and segment combination.

    Returns
    -------
    None
    """
    
    segment, hvlevel = args['segment'], args['hv_lvl']

    settings = get_settings()
    database = get_database()
    database.connect()
    
    #-- Nothing bad before 2010,
    #-- and there are some weird gainmaps back there
    #-- filtering out for now (expstart > 55197).
    
    #-- Gain.y < 300 and Gain.y > 200 are counts in the spectral region.
    #-- Don't care about anything outside of that. 
    all_coords = list(Gain.select().distinct().where(
                                                (Gain.hv_lvl == hvlevel) &
                                                (Gain.gain <= 3) &
                                                (Gain.counts >= 30) &
                                                (Gain.segment == segment) &
                                                (Gain.expstart > 55197) &
                                                (Gain.y > 400//Y_BINNING) &
                                                (Gain.y < 600//Y_BINNING) &
                                                Gain.filename.not_in(Flagged_Pixels.select(Flagged_Pixels.filename))
                                                ).dicts())
    
    #-- Set up collection dictionary.
    result = collections.defaultdict(list)
    
    #-- for each x,y coodinate pair, create a dictionary where there key is the x,y pair that
    #-- contains a list of dictionaries for all of the entries of that x,y pair.
    for d in all_coords:
        result[d['x'], d['y']].append(d)
    
    #-- For each key in the dictionary, find the dictionary that has the date entry where the pixel
    #-- first went bad.
    bad_pix = []
    for coord_pair in result.keys():
        row_bad = min(result[coord_pair], key=lambda x:x['expstart'])
        
        #-- Take meta and make dict according to fields in Flagged_Pixels and append to list.
        bad_dict = {'segment': segment,
                    'hv_lvl': hvlevel,
                    'x': row_bad['x'],
                    'y': row_bad['y'],
                    'filename': row_bad['filename'],
                    'mjd_below_3': row_bad['expstart']}
        
        bad_pix.append(bad_dict)

    bulk_insert(Flagged_Pixels, bad_pix)
    database.close()
#-------------------------------------------------------------------------------
def populate_bad_pix():
    """
    Main driver to populate flagged pixel table. 

    Parameters
    ----------
    None

    Retruns
    -------
    None
    """

    #-- Connect to DB
    settings = get_settings()
    database = get_database()
    database.connect()

    #-- Get all of the possible segments and HV level combinations
    all_combos = Gain.select(Gain.segment, Gain.hv_lvl).distinct().dicts()
    database.close()

    #-- Populate the database.
    pool = mp.Pool(processes=settings['num_cpu'])
    pool.map(find_flagged, all_combos)
#-------------------------------------------------------------------------------