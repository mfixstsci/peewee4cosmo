
from __future__ import absolute_import

"""Routine to monitor the modal gain in each pixel as a
function of time.  Uses COS Cumulative Image (CCI) files
to produce a modal gain map for each time period.  Modal gain
maps for each period are collated to monitor the progress of
each pixel(superpixel) with time.  Pixels that drop below
a threshold value are flagged and collected into a
gain sag table reference file (gsagtab).

The PHA modal gain threshold is set by global variable MODAL_GAIN_LIMIT.
Allowing the modal gain of a distribution to come within 1 gain bin
of the threshold results in ~8% loss of flux.  Within
2 gain bins, ~4%
3 gain bins, ~2%
4 gain bins, ~1%

However, due to the column summing, a 4% loss in a region does not appear to be so in the extracted spectrum.
"""

__author__ = 'Justin Ely, Mees Fix'
__maintainer__ = 'Mees Fix'
__email__ = 'mfix@stsci.edu'
__status__ = 'Active'

import os
import shutil
import time
from datetime import datetime
import glob
import sys
from sqlalchemy.engine import create_engine
import logging
logger = logging.getLogger(__name__)

from astropy.io import fits
from astropy.time import Time

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sqlalchemy.sql.functions import concat

from ..utils import send_email
from .constants import *  

from ..database.models import get_database, get_settings, Files
from ..database.models import Flagged_Pixels, Gain
from .gainmap import make_all_gainmaps
from .gainmap_sagged_pixel_overplotter import make_overplot, hotspot_plotter_interactive, gsagtab_overplot_comparison

import collections
import functools
import multiprocessing as mp
#------------------------------------------------------------

def main(out_dir, hotspot_filter=True):
    """Main driver for monitoring program.

    Parameters
    ----------
    out_dir: str
        Strings where you would like to write files to
    
    Returns
    -------
    None
    """
    
    settings = get_settings()

    logger.info("STARTING MONITOR")
    
    logger.info("MAKING HOTSPOT PLOTS")
    hotspot_plotter_interactive('FUVA')
    hotspot_plotter_interactive('FUVB')
    
    logger.info("MAKING NEW GSAGTAB")
    reg_gsagtab = make_gsagtab_db(out_dir, filter=hotspot_filter)
    blu_gsagtab = make_gsagtab_db(out_dir, blue=True)
    gsagtabs = [reg_gsagtab, blu_gsagtab]

    #-- Pooled gainmap creation.
    logger.info("MAKING COMBINED GAINMAPS")
    partial = functools.partial(make_all_gainmaps, 
                                gainmap_dir=os.path.join(settings['monitor_location'],'CCI'))
    
    pool = mp.Pool(processes=settings['num_cpu'])
    pool.map(partial, range(150,179))

    #-- HV Level doest matter when total=True.    
    make_all_gainmaps(100, gainmap_dir=os.path.join(settings['monitor_location'],'CCI'), start_mjd=55055, end_mjd=70000, total=True)
    
    # logger.info("MAKING GAINMAP + GSAG OVERPLOT")
    # make_overplot(reg_gsagtab)
    
    #-- Current CRDS gsagtabs
    current_gsag_tab = os.path.join(settings['lref'], 'zbn1927gl_gsag.fits')    
    current_blue_tab = os.path.join(settings['lref'], 'zbn1927fl_gsag.fits')
    
    cdbs_gsagtabs = [current_gsag_tab, current_blue_tab]
    #-- Create comparison figures for HV 163-175
    pool = mp.Pool(processes=settings['num_cpu'])
    
    for new_gsagtab, calcos_tab in zip(gsagtabs, cdbs_gsagtabs):
        #-- Set up partials for multiprocessing.
        partial = functools.partial(gsagtab_overplot_comparison,
                                    compare=True, 
                                    potential_gsagtab=new_gsagtab,
                                    current_gsagtab=calcos_tab)
        
        pool.map(partial, [163,167,169,171,173,175,178])
        
        #-- Just make overplots without comparing.
        partial = functools.partial(gsagtab_overplot_comparison,
                                    potential_gsagtab=new_gsagtab)
        pool.map(partial, [163,167,169,171,173,175,178])

    logger.info("FINISH MONITOR")

#------------------------------------------------------------

def gsagtab_extension(date, lx, dx, ly, dy, dq, dethv, hv_string, segment):
    """Creates a properly formatted gsagtab table from input columns.

    Parameters
    ----------
    date: list
        List of dates pixels went bad.
    lx: list
        List of x pixel positions (xcorr) space.
    dx: list
        List of x bin for COS CCI super pixels.
    ly: list
        List of y pixel positions (ycorr) space.
    dy: list
        List of y bin for COS CCI super pixels.
    dq: list
        List of dq values (8192) for each flagged pixel
    dethv: list
        List of voltage values for flagged pixels.
    hv_string: str
        HVLEVELA or HLEVELB based on FUVA or FUVB for different header exts.
    segment: str
        FUVA or FUVB
    """

    #-- Set arrays
    lx = np.array(lx)
    ly = np.array(ly)
    dx = np.array(dx)
    dy = np.array(dy)
    dq = np.array(dq)

    #-- Create data columns for fits file.
    date_col = fits.Column('DATE','D','MJD',array=date)
    lx_col = fits.Column('LX','J','pixel',array=lx)
    dx_col = fits.Column('DX','J','pixel',array=dx)
    ly_col = fits.Column('LY','J','pixel',array=ly)
    dy_col = fits.Column('DY','J','pixel',array=dy)
    dq_col = fits.Column('DQ','J','',array=dq)
    
    #-- Create data table
    tab = fits.TableHDU.from_columns([date_col,lx_col,ly_col,dx_col,dy_col,dq_col])

    #-- Add comments.
    tab.header.add_comment(' ',after='TFIELDS')
    tab.header.add_comment('  *** Column formats ***',after='TFIELDS')
    tab.header.add_comment(' ',after='TFIELDS')
    tab.header.set(hv_string, dethv, after='TFIELDS',comment='High voltage level')
    tab.header.set('SEGMENT', segment, after='TFIELDS')
    tab.header.add_comment(' ',after='TFIELDS')
    tab.header.add_comment('  *** End of mandatory fields ***',after='TFIELDS')
    tab.header.add_comment(' ',after='TFIELDS')

    return tab

#------------------------------------------------------------

def date_string(date_time):
    """Takes a datetime object and returns
    a pedigree formatted string.
    
    Parameters
    ----------
    date_time: str
        Date when monitor is run.

    Returns
    -------
    date_string: str
        A formatted string used in filenames and header keywords.
    
    """

    day = str(date_time.day)
    month = str(date_time.month)
    year = str(date_time.year)

    if len(day) < 2:
        day = '0' + day

    if len(month) < 2:
        month = '0' + month

    date_string = day + '/' + month + '/' + year
    
    return date_string
#------------------------------------------------------------

def in_boundary(segment, ly, dy):
    """
    Check to see if LP2 blue mode pixel is in LP1 gsag region
    and exclude it from the the blue mode GSAGTAB.

    Parameters
    ----------
    segement: str
        FUVA or FUVB
    ly: int
        y pixel position
    dy: int
        delta y for binning.

    Returns
    -------
    True or False
    """

    #-- Roughly the area around
    boundary = {'FUVA': 493, 'FUVB': 557}
    padding = 4

    boundary_pix = set(np.arange(boundary[segment]-padding,
                                 boundary[segment]+padding+1))

    affected_pix = set(np.arange(ly, ly+dy+1))

    if affected_pix.intersection(boundary_pix):
        return True

    return False

#------------------------------------------------------------

def make_gsagtab_db(out_dir, blue=False, filter=False):
    """Creates GSAGTAB from Flagged_Pixel DB table.

    Parameters
    ----------
    out_dir: str
        Directory to write figures/files out to
    blue: bool
        If true make bluemode maps
    filter: bool
        If filter, check pixels for temporary sagging
   
    Returns
    -------
    None

    Products
    --------
    new_gsagtab.fits
    """

    if filter:
        logger.info('CHECKING LP PIXEL RECOVERY')
        check_pixel_recovery('FUVA')
        check_pixel_recovery('FUVB')
        filename = 'gsag_filter_%s.fits'%(TIMESTAMP)
    else:
        filename = 'gsag_%s.fits'%(TIMESTAMP)
    
    out_fits = os.path.join(out_dir, filename)

    #-- Begin header data.
    hdu_out=fits.HDUList(fits.PrimaryHDU())
    date_time = str(datetime.now())
    date_time = date_time.split()[0]+'T'+date_time.split()[1]
    hdu_out[0].header['DATE'] = (date_time, 'Creation UTC (CCCC-MM-DD) date')
    hdu_out[0].header['TELESCOP'] = 'HST'
    hdu_out[0].header['INSTRUME'] = 'COS'
    hdu_out[0].header['DETECTOR'] = 'FUV'
    hdu_out[0].header['COSCOORD'] = 'USER'
    hdu_out[0].header['VCALCOS'] = '2.0'
    hdu_out[0].header['USEAFTER'] = 'May 11 2009 00:00:00'
    hdu_out[0].header['CENWAVE'] = 'N/A'

    today_string = date_string(datetime.now())
    hdu_out[0].header['PEDIGREE'] = 'INFLIGHT 25/05/2009 %s'%(today_string)
    hdu_out[0].header['FILETYPE'] = 'GAIN SAG REFERENCE TABLE'

    descrip_string = 'Gives locations of gain-sag regions as of %s'%( str(datetime.now().date() ))
    while len(descrip_string) < 67:
        descrip_string += '-'
    hdu_out[0].header['DESCRIP'] = descrip_string
    hdu_out[0].header['COMMENT'] = ("= 'This file was created by M. Fix'")
    hdu_out[0].header.add_history('Flagged regions in higher voltages have been backwards populated')
    hdu_out[0].header.add_history('to all lower HV levels for the same segment.')
    hdu_out[0].header.add_history('')
    hdu_out[0].header.add_history('A region will be flagged as bad when the detected')
    hdu_out[0].header.add_history('flux is found to drop by 5%.  This happens when')
    hdu_out[0].header.add_history('the measured modal gain of a region falls to ')
    hdu_out[0].header.add_history('%d given current lower pulse height filtering.'%(MODAL_GAIN_LIMIT) )

    #-- Possible HVs for all of the extensions.
    possible_hv_strings = ['000', '100'] + list(map(str, list(range(142, 179))))

    #-- Connect to database.
    settings = get_settings()
    database = get_database()
    database.connect()

    #-- Get all segments
    segment_results = Gain.select(Gain.segment).distinct()

    #-- Put segments in list.
    segments = [row.segment for row in segment_results]

    #-- For each segment, loop through all the possible HVs and create an extension in the GSAGTAB.
    for segment in segments:
        hvlevel_string = 'HVLEVEL' + segment[-1].upper()

        #-- Begin HV loop.
        for hv_level in possible_hv_strings:
            #-- Set all of the empty lists that will be the data arrays for GSAGTAB exts.
            date = []
            lx = []
            dx = []
            ly = []
            dy = []
            dq = []

            hv_level = int(hv_level)

            #-- Logic to filter hotspots or nah.
            if filter:
                #-- Get all of the sagged pixels for a specific Segment/HV Level combo.
                flagged_pix = list(Flagged_Pixels.select().distinct().where(
                                                                            (Flagged_Pixels.segment == segment) &
                                                                            (Flagged_Pixels.hv_lvl >= hv_level) &
                                                                            (Flagged_Pixels.recovery == False)
                                                                            ).dicts())
            else:
                flagged_pix = list(Flagged_Pixels.select().distinct().where(
                                                                            (Flagged_Pixels.segment == segment) &
                                                                            (Flagged_Pixels.hv_lvl >= hv_level)
                                                                            ).dicts())
            result = collections.defaultdict(list)
            #-- for each x,y coodinate pair, create a dictionary where there key is the x,y pair that
            #-- contains a list of dictionaries for all of the entries of that x,y pair.
            for d in flagged_pix:
                result[d['x'], d['y']].append(d)
            
            #-- Find the min date for each coordinate pair.
            bad_pix = []
            for coord_pair in result.keys():
                row_bad = min(result[coord_pair], key=lambda x:x['mjd_below_3'])

                #-- Block for blue mode creation. in_boundary checks for LP1 gsag holes that are in the
                #-- LP2 extraction regions. We want to exclude them because we get low counts in the wings
                #-- and since blue modes use boxcar extraction it will flag the whole column and we want to avoid that.
                #-- According to the team, the long term solution would be to only have one gsagtab and LP2 will use 2zone extraction.
                if blue and in_boundary(row_bad['segment'], row_bad['y']*Y_BINNING, Y_BINNING):
                    logger.debug("EXCLUDING FOR BLUE MODES: {} {} {}".format(row_bad['segment'], row_bad['y']*Y_BINNING, Y_BINNING))
                    continue

                #-- Apply binning constant and append values to lists
                lx.append(row_bad['x']*X_BINNING)
                dx.append(X_BINNING)
                ly.append(row_bad['y']*Y_BINNING)
                dy.append(Y_BINNING)
                date.append(row_bad['mjd_below_3'])
                dq.append(8192)

            if not len(lx):
                #-- Extension tables cannot have 0 entries, a
                #-- region of 0 extent centered on (0,0) is
                #-- sufficient to prevent CalCOS crash.
                lx.append(0)
                ly.append(0)
                dx.append(0)
                dy.append(0)
                date.append(0)
                dq.append(8192)

            logger.debug('found {} bad regions'.format(len(date)))
            tab = gsagtab_extension(date, lx, dx, ly, dy, dq, hv_level, hvlevel_string, segment)
            hdu_out.append(tab)
    
    #-- If blue mode, change the naming scheme and the header.
    if blue:
        out_fits = out_fits.replace('.fits', '_blue.fits')
        hdu_out[0].header['CENWAVE'] = 'BETWEEN 1055 1097'

        descrip_string = 'Blue-mode gain-sag regions as of %s'%(str(datetime.now().date()))
        while len(descrip_string) < 67:
            descrip_string += '-'
        hdu_out[0].header['DESCRIP'] = descrip_string

    database.close()
    hdu_out.writeto(out_fits, clobber=True)
    logger.info('WROTE: GSAGTAB to %s'%(out_fits))
    return out_fits

#------------------------------------------------------------

def check_pixel_recovery(segment):
    """Check pixels for actual sagging. Sometimes
    a pixel can drop below gain 3 but can recover.

    Parameters
    ----------
    segment: str
        FUVA or FUVB
    
    Returns
    -------
    None
    """

    #-- Connect to database
    settings = get_settings()
    database = get_database()
    database.connect()
    
    #-- Profile widths for LP4. (Super Pixel binned by two.)
    lp4_profile ={'FUVA':[400//Y_BINNING, 440//Y_BINNING],
                  'FUVB':[460//Y_BINNING, 500//Y_BINNING]}


    #-- Get all of the flagged pixels at LP4
    flagged_pixels = Flagged_Pixels.select().distinct().where(
                                                             (Flagged_Pixels.recovery == False) &
                                                             (Flagged_Pixels.y.between(lp4_profile[segment][0],lp4_profile[segment][1])) &
                                                             (Flagged_Pixels.x.between(1000//X_BINNING,15000//X_BINNING)) &
                                                             (Flagged_Pixels.segment == segment)
                                                             )
    
    #-- Get all measurements for all of the x,y pixels returned.    
    gain_measurements = list(Gain.select().distinct().where(
                                                           (Gain.x.in_([row.x for row in flagged_pixels])) &
                                                           (Gain.y.in_([row.y for row in flagged_pixels])) &
                                                           (Gain.segment == segment)
                                                           ).dicts())
    
    logger.info('FOUND {} SAGGED PIXELS FOR SEGMENT {}!'.format(len(gain_measurements),segment))
    
    result = collections.defaultdict(list)
    
    #-- Organize all of the pixels by segment.
    for d in gain_measurements:
        result[d['segment'], d['x'], d['y']].append(d)


    #-- For every segment,x,y combo check the latest measurements
    for combo in result.keys():
        segment,x,y = combo
        
        #-- Get the latest measurement for this pixel.
        latest_measurement = max(result[combo], key=lambda x:x['expstart'])
        
        #-- If the latest gain measurement is above 3, the pixel recovered!
        if latest_measurement['gain'] > 3:
            Flagged_Pixels.update(recovery=True).where(
                                                      (Flagged_Pixels.x == latest_measurement['x']) &
                                                      (Flagged_Pixels.y == latest_measurement['y']) &
                                                      (Flagged_Pixels.segment == latest_measurement['segment'])
                                                      ).execute()
            
            logger.info('UPDATING ROW FOR x={}, y={}! gain={} AS OF {}'.format(latest_measurement['x']*X_BINNING,
                                                                               latest_measurement['y']*Y_BINNING,
                                                                               latest_measurement['gain'],
                                                                               Time(latest_measurement['expstart'],format='mjd').iso)) 
        #-- Else, the pixel never recovered keep the entry.
        else:
            logger.warning('PIXEL x={}, y={} IS SAGGED! gain={} AS OF {}.'.format(latest_measurement['x']*X_BINNING,
                                                                                  latest_measurement['y']*Y_BINNING,
                                                                                  latest_measurement['gain'],
                                                                                  Time(latest_measurement['expstart'],format='mjd').iso))
            continue

    database.close()