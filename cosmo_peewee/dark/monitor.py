"""Perform regular monitoring of the COS FUV and NUV dark rates
"""

from __future__ import print_function, absolute_import, division

import os
import datetime
import numpy as np
import shutil
import glob
import math
import logging
logger = logging.getLogger(__name__)

from astropy.io import fits
from astropy.time import Time

from calcos import orbit
from calcos.timeline import gmst, ASECtoRAD, DEGtoRAD, eqSun, DIST_SUN, RADIUS_EARTH, computeAlt, computeZD, rectToSph

from .solar import get_solar_data
from .plotting import plot_histogram, plot_time, plot_orbital_rate
from .interactive_plots import plot_time as interactive_plot_time

from ..utils import corrtag_image, remove_if_there
from ..database.models import get_settings, get_database
from ..database.models import Darks

from copy import deepcopy

#-------------------------------------------------------------------------------

def get_sun_loc(mjd, full_path):
    
    """ Get the location of the sun from SPT files. 

    Parameters
    ----------
    mjd : float
        The MJD for an exposure
    full_path : str
        String of path and filename

    Yields
    ------
    long_sun : float
        Longitude of the sun
    lat_sun : float
        Latitude of the sun
    """

    rootname = fits.getval(full_path, 'ROOTNAME')

    path, _ = os.path.split(full_path)
    sptfile = os.path.join(path, rootname + '_spt.fits')
    if not os.path.exists(sptfile):
        sptfile += '.gz'
        if not os.path.exists(sptfile):
            raise IOError("Cannot find sptfile {}".format(sptfile))

    orb = orbit.HSTOrbit(sptfile)

    if isinstance(mjd, (int, float)):
        mjd = list(mjd)

    for m in mjd:
        (rect_hst, vel_hst) = orb.getPos(m)
        (r, ra_hst, dec_hst) = rectToSph(rect_hst)

        #-- Assume that we want geocentric latitude.  The difference from
        #-- astronomical latitude can be up to about 8.6 arcmin.
        lat_hst = dec_hst
        
        #-- Subtract the sidereal time at Greenwich to convert to longitude.
        long_hst = ra_hst - 2. * math.pi * gmst(m)
        if long_hst < 0.:
            long_hst += (2. * math.pi)

        long_col = long_hst / DEGtoRAD
        lat_col = lat_hst / DEGtoRAD
        
        #-- equatorial coords of the Sun
        rect_sun = eqSun(m)

        (r, ra_sun, dec_sun) = rectToSph(rect_sun)
        lat_sun = dec_sun
        long_sun = ra_sun - 2. * math.pi * gmst(m)
        if long_sun < 0.:
            long_sun += (2. * math.pi)
        long_sun /= DEGtoRAD
        lat_sun /= DEGtoRAD

        yield long_sun, lat_sun

#-------------------------------------------------------------------------------

def get_temp(filename):
    
    """Get detector temperture during observation from spt filename

    Parameters
    ----------
    filename : str
        FITS file for which the temperature is to be Found

    Returns
    -------
    temperature : float
        Detector temperature at the time of the observation

    """

    with fits.open(filename) as hdu:
        detector = hdu[0].header['DETECTOR']
        segment = hdu[0].header['SEGMENT']
        rootname = hdu[0].header['ROOTNAME']

    if detector == 'FUV' and segment == 'FUVA':
        temp_keyword = 'LDCAMPAT'
    elif detector == 'FUV' and segment == 'FUVB':
        temp_keyword = 'LDCAMPBT'
    elif detector == 'NUV':
        temp_keyword = 'LMMCETMP'
    else:
        raise ValueError('WHAT DETECTOR AND SEGMENTS ARE THESE?! {} {}'.format(detector, segment))

    path, name = os.path.split(filename)
    spt_file = os.path.join(path, rootname + '_spt.fits')

    try:
        temperature = fits.getval(spt_file, temp_keyword, ext=2)
    except IOError:
        temperature = fits.getval(spt_file + '.gz', temp_keyword, ext=2)

    return temperature

#-------------------------------------------------------------------------------

def mjd_to_decyear(time_array):
    
    """ Changes the date in MJD units to decimal years.
    
    Parameters
    ----------
    time_array : array like
        A list of times measured in MJD

    Returns
    -------
    out_times : np.array
        A numpy array of MJD to decimal year conversions.
    """

    times = Time(time_array, scale='tt', format='mjd')

    out_times = []
    for value in times:
        year = value.datetime.year
        n_days = (value.datetime - datetime.datetime(value.datetime.year, 1, 1)).total_seconds()
        total_days = (datetime.datetime(value.datetime.year+1, 1, 1) - datetime.datetime(value.datetime.year, 1, 1)).total_seconds()

        fraction = float(n_days) / total_days

        out_times.append(year + fraction)

    return np.array(out_times)

#-------------------------------------------------------------------------------

def pull_orbital_info(data_object, step=25):
    
    """ Pull second by second orbital information. This function populates the
    darks DB table.

    Parameters
    ----------
    data_object : peewee query result
        Contains the path and filename information needed to process data.
    step : int
        Time step in seconds

    Yields
    ------
    info : dictionary
        A dictionary of meta data that will be stored in a row of the DB table. 

    """
    
    full_path = os.path.join(data_object.path, data_object.filename)
    
    SECOND_PER_MJD = 1.15741e-5

    info = {}
    info['filename'] = os.path.split(full_path)[-1]

    hdu = fits.open(full_path)
    
    #-- Get timeline extension from corrtag
    try:
        timeline = hdu['timeline'].data
        segment = hdu[0].header['segment']
    except KeyError:
        logger.debug("NO TIMELINE EXTENSION FOUND FOR: {}".format(full_path))
        print(info)
        yield info
        raise StopIteration

    #-- Set boundaries based on segment/detector
    if segment == 'N/A':
        segment = 'NUV'
        xlim = (0, 1024)
        ylim = (0, 1204)
        pha = (-1, 1)
    elif segment == 'FUVA':
        xlim = (1200, 15099)
        ylim = (380, 680)
        pha = (2, 23)
    elif segment == 'FUVB':
        xlim = (950, 15049)
        ylim = (440, 720)
        pha = (2, 23)
    else:
        raise ValueError('WHAT SEGMENT IS THIS?! {}'.format(segment))

    #-- Get some basic header info.
    info['rootname'] = hdu[0].header['rootname']
    info['targname'] = hdu[0].header['targname']
    info['detector'] = segment
    info['temp'] = get_temp(full_path)

    #-- Break timeline into 25 second intervals from exstart to end.
    times = timeline['time'][::step].copy()

    #-- Get all of the latitude/longtitude measurements.
    lat = timeline['latitude'][:-1][::step].copy().astype(np.float64)
    lon = timeline['longitude'][:-1][::step].copy().astype(np.float64)

    #-- Convert time from mjd.
    mjd_per_step = hdu[1].header['EXPSTART'] + \
                   times.copy().astype(np.float64) * \
                   SECOND_PER_MJD
    
    #-- Get the solar longitude and latitude.
    sun_lat = []
    sun_lon = []
    for item in get_sun_loc(mjd_per_step, full_path):
        sun_lon.append(item[0])
        sun_lat.append(item[1])

    #-- Grab the last value in the list, this value is used in the plot to show the variation
    #-- during the exposure.
    mjd = mjd_per_step[:-1]

    #-- Convert from mjd to decimal year.
    decyear = mjd_to_decyear(mjd)

    #-- Make sure that the exposure actually happened.
    if not len(times):
        logger.debug("TIME ARRAY EMPTY FOR: {}".format(full_path))
        blank = np.array([0])
        print(info)
        yield info
        raise StopIteration

    #-- Grab the events list and filter based on PHA that is kept.
    events = hdu['events'].data
    filtered_index = np.where((events['PHA'] > pha[0]) &
                              (events['PHA'] < pha[1]) &
                              (events['XCORR'] > xlim[0]) &
                              (events['XCORR'] < xlim[1]) &
                              (events['YCORR'] > ylim[0]) &
                              (events['YCORR'] < ylim[1]))

    #-- For TA's we keep all of the PHA extensions.
    ta_index = np.where((events['XCORR'] > xlim[0]) &
                        (events['XCORR'] < xlim[1]) &
                        (events['YCORR'] > ylim[0]) &
                        (events['YCORR'] < ylim[1]))

    #-- Build a histogram of the events based on PHA extensions.
    counts = np.histogram(events[filtered_index]['time'], bins=times)[0]
    ta_counts = np.histogram(events[ta_index]['time'], bins=times)[0]

    #-- Calculate dark rate [counts/pix/sec]
    npix = float((xlim[1] - xlim[0]) * (ylim[1] - ylim[0]))
    counts = counts / npix / step
    ta_counts = ta_counts / npix / step

    #-- Make sure arrays have same length so we dont crash.
    if not len(lat) == len(counts):
        lat = lat[:-1]
        lon = lon[:-1]
        sun_lat = sun_lat[:-1]
        sun_lon = sun_lon[:-1]

    assert len(lat) == len(counts), \
        'ARRAYS ARE NOT EQUAL LENGTH {}:{}'.format(len(lat), len(counts))

    #-- Yield results for each time step that get added to DB.
    if not len(counts):
        logger.debug("ZERO-LENGTH ARRAY FOUND FOR: {}".format(full_path))
        yield info
    else:
        for i in range(len(counts)):
            #-- better solution than round?
            info['date'] = round(decyear[i], 3)
            info['dark'] = round(counts[i], 7)
            info['ta_dark'] = round(ta_counts[i], 7)
            info['latitude'] = round(lat[i], 7)
            info['longitude'] = round(lon[i], 7)
            info['sun_lat'] = round(sun_lat[i], 7)
            info['sun_lon'] = round(sun_lon[i], 7)
            info['mjd_per_step'] = mjd_per_step[i]
            
            yield deepcopy(info)

#-------------------------------------------------------------------------------

def compile_phd():

    #-- THIS STILL USES SQLALCHEMY FROM cos_monitoring.
    #-- MAY DELETE IN THE FUTURE.
    raise NotImplementedError("Nope, seriously can't do any of this.")

    #-- populate PHD table

    columns = ', '.join(['bin{} real'.format(pha) for pha in range(0,31)])
    c.execute("""CREATE TABLE {} ( obsname text, {})""".format(table, columns ))


    c.execute( """SELECT obsname FROM %s """ %(table))
    already_done = set( [str(item[0]) for item in c] )

    for filename in available:
        obsname = os.path.split(filename)[-1]
        if obsname in already_done:
            print(filename, 'done')
        else:
            print(filename, 'running')

        counts = pha_hist(filename)
        table_values = (obsname, ) + tuple(list(counts) )

        c.execute( """INSERT INTO %s VALUES (?{})""" % (table, ',?'*31 ),
                   table_values)

        db.commit()

#-------------------------------------------------------------------------------

def pha_hist(filename):
    hdu = fits.open( filename )
    pha_list_all = hdu[1].data['PHA']
    counts, bins = np.histogram(pha_list_all, bins=31, range=(0, 31))

    return counts

#-------------------------------------------------------------------------------

def make_plots(detector, base_dir, TA=False, mjd_per_step=False):
    
    """ Create static monitoring plots for FUV/NUV dark rates.

    Parameters
    ----------
    detector : str
        The COS mode trends you are interested in plotting.
    base_dir : str
        Directory you are interested in writing to.
    TA : bool
        Flag to monitor target acq dark rate.

    Returns
    -------
    None
    """

    if detector == 'FUV':
        search_strings = ['_corrtag_a.fits', '_corrtag_b.fits']
        segments = ['FUVA', 'FUVB']
    elif detector == 'NUV':
        search_strings = ['_corrtag.fits']
        segments = ['NUV']
    else:
        raise ValueError('Only FUV or NUV allowed.  NOT:{}'.format(detector) )

    try:
        solar_data = np.genfromtxt(os.path.join(base_dir, 'solar_flux.txt'), dtype=None)
        if mjd_per_step:
            solar_date = np.array([float(str(line[0]).replace(" ", "-")) for line in solar_data])
        else:
            solar_date = np.array(mjd_to_decyear([line[0] for line in solar_data]))
        solar_flux = np.array([line[1] for line in solar_data])
    except TypeError:
        logger.warning("COULDN'T READ SOLAR DATA. PUTTING IN ZEROS.")
        solar_date = np.ones(1000)
        solar_flux = np.ones(1000)

    #-- Open settings and get database
    settings = get_settings()
    database = get_database()


    if TA:
        logger.info("MAKING PLOTS FOR {} TARGACQS".format(detector))
    else:
        logger.info("MAKING PLOTS FOR {}".format(detector))

    for key, segment in zip(search_strings, segments):
        
        logger.debug('CREATING TIME PLOT FOR {}:{}'.format(segment, key))
        #-- Query for data here!
        data = Darks.select().where(Darks.detector == segment)

        #-- Parse whether you want to plot dark monitoring or targacq dark.
        if TA:
            dark_key = 'ta_dark'
            dark = np.array([item.ta_dark for item in data])
        else:
            dark_key = 'dark'
            dark = np.array([item.dark for item in data])
        
        #-- Break query into all of it's components.
        temp = np.array([item.temp for item in data])
        latitude = np.array([item.latitude for item in data])
        longitude = np.array([item.longitude for item in data])
        sun_lat = np.array([item.sun_lat for item in data])
        sun_lon = np.array([item.sun_lon for item in data])
        if mjd_per_step:
            mjd = np.array([item.mjd_per_step for item in data])
        else:
            mjd = np.array([item.date for item in data])

        #-- Sort data by date.
        index = np.argsort(mjd)
        mjd = mjd[index]
        dark = dark[index]
        temp = temp[index]
        latitude = latitude[index]
        longitude = longitude[index]
        sun_lat = sun_lat[index]
        sun_lon = sun_lon[index]
        date = np.array([item.date for item in data])

        #-- Plot vs orbit
        logger.debug('CREATING ORBIT PLOT FOR {}:{}'.format(segment, key))
        
        outname = os.path.join(base_dir, detector, '{}_vs_orbit_{}.png'.format(dark_key, segment))
        plot_orbital_rate(longitude, latitude, dark, sun_lon, sun_lat, outname)

        #-- Plot histogram of darkrates
        logger.debug('CREATING HISTOGRAM PLOT FOR {}:{}'.format(segment, key))

        if not mjd_per_step:
            for year in set(map(int, date)):
                if year == 0:
                    continue
                else:
                    index = np.where( (date >= year) &
                                    (date < year + 1))
                    hist_outname = os.path.join(base_dir, detector, '{}_hist_{}_{}.png'.format(dark_key, year, segment))
                    plot_histogram(dark[index], hist_outname)

                    inter_outname = os.path.join(base_dir, detector, '{}_interactive_{}_{}.html'.format(dark_key, year, segment))
                    interactive_plot_time(detector, dark[index], mjd[index], temp[index], solar_flux, solar_date, inter_outname)

            index = np.where(date >= date.max() - .5)
            outname = os.path.join(base_dir, detector, '{}_hist_-6mo_{}.png'.format(dark_key, segment))
            plot_histogram(dark[index], outname )

            outname = os.path.join(base_dir, detector, '{}_hist_{}.png'.format(dark_key, segment))
            plot_histogram(dark, outname)

        #-- Dark vs Time Plots.
        #-- Restrict to avoid SAA (?)
        index_keep = np.where((longitude < 250) | (latitude > 10))[0]
        mjd = mjd[index_keep]
        dark = dark[index_keep]
        temp = temp[index_keep]
        
        logger.debug('CREATING INTERACTIVE VS TIME PLOT FOR {}:{}'.format(segment, key))
        
        #-- Interactive plots
        if mjd_per_step:
            outname = os.path.join(base_dir, detector, '{}_vs_step_time_{}.html'.format(dark_key, segment))
        else:
            outname = os.path.join(base_dir, detector, '{}_vs_time_{}.html'.format(dark_key, segment))
        
        interactive_plot_time(detector, dark, mjd, temp, solar_flux, solar_date, outname)

        if not mjd_per_step:
            #-- Static plots
            if mjd_per_step:
                outname = os.path.join(base_dir, detector, '{}_vs_step_time_{}.png'.format(dark_key, segment))
            else:
                outname = os.path.join(base_dir, detector, '{}_vs_time_{}.png'.format(dark_key, segment))

            plot_time(detector, dark, mjd, temp, solar_flux, solar_date, outname)
#-------------------------------------------------------------------------------

def move_products(base_dir, web_dir):
    
    '''Move monitoring figures to webpage directory. 
    
    Parameters
    ----------
    base_dir : str
        Directory where figures are located
    web_dir : 
        COS monitoring webpage directory.

    Returns
    -------
    None
    '''

    for detector in ['FUV', 'NUV']:

        #-- Where would you like to write the plots to?
        write_dir = os.path.join(web_dir, detector.lower() + '_darks/')
        
        #-- If the path doesnt exist, make your own...
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
        
        #-- Combine the base monitoring dir with the detector specific dir.
        detector_dir = os.path.join(base_dir, detector)
        
        #-- Grab all of the files you wish to move....
        move_list = glob.glob(detector_dir + '/*.p??')
                
        for item in move_list:
            try:
                #-- Don't want any python scripts moving.
                if item.endswith('.py~'):
                    logger.debug("REMOVING {}".format(item))
                    move_list.remove(item)
                    continue
                else:
                    logger.debug("MOVING {}".format(item))
                
                #-- Split the file and paths.
                path, file_to_move = os.path.split(item)
                
                #-- Update the permissions
                os.chmod(item, 0o766)

                #-- Remove the file if it exists in the webpage dir.
                remove_if_there(write_dir + file_to_move)
                
                #-- Copy the file over.
                shutil.copy(item, write_dir + file_to_move)

            except OSError:
                logger.warning("HIT AN OS ERROR FOR {}, LEAVING IT THERE".format(item))
                move_list.remove(item)
        
        #-- Change all of the file permissions.
        os.system('chmod 777 ' + write_dir + '/*.png')
#-------------------------------------------------------------------------------

def monitor():
    """Main monitoring pipeline
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    
    """

    logger.info("STARTING MONITOR")

    settings = get_settings()
    out_dir = os.path.join(settings['monitor_location'], 'Darks')
    web_dir = settings['webpage_location']

    if not os.path.exists(out_dir):
        logger.warning("CREATING OUTPUT DIRECTORY: {}".format(out_dir))
        os.makedirs(out_dir)

    get_solar_data(out_dir)

    for detector in ['FUV', 'NUV']:
        make_plots(detector, out_dir)
        make_plots(detector, out_dir, mjd_per_step=True)
        if detector == 'FUV':
            make_plots(detector, out_dir, TA=True)
            make_plots(detector, out_dir, TA=True, mjd_per_step=True)

    logger.info("MOVING PRODUCTS TO WEB DIRECTORY")
    move_products(out_dir, web_dir)

#-------------------------------------------------------------------------------
