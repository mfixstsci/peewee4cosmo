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

from ..utils import corrtag_image
from ..database.models import get_settings, get_database
from ..database.models import Darks

from sqlalchemy.sql.functions import concat

#-------------------------------------------------------------------------------

def get_sun_loc(mjd, full_path):
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
        raise ValueError('What??? {} {}'.format(detector, segment))

    path, name = os.path.split(filename)
    spt_file = os.path.join(path, rootname + '_spt.fits')

    try:
        temperature = fits.getval(spt_file, temp_keyword, ext=2)
    except IOError:
        temperature = fits.getval(spt_file + '.gz', temp_keyword, ext=2)

    return temperature

#-------------------------------------------------------------------------------

def mjd_to_decyear( time_array ):
    """ pull this out when you get it into astropy.time

    """

    times = Time( time_array, scale='tt', format='mjd' )

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
    """ Pull second by second orbital information from the full_path

    """
    
    full_path = os.path.join(data_object.path, data_object.filename)
    
    SECOND_PER_MJD = 1.15741e-5

    info = {}
    info['filename'] = os.path.split(full_path)[-1]

    hdu = fits.open(full_path)
    
    try:
        timeline = hdu['timeline'].data
        segment = hdu[0].header['segment']
    except KeyError:
        logger.debug("no timeline extension found for: {}".format(full_path))
        yield info
        raise StopIteration

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
        raise ValueError('What segment is this? {}'.format(segment))

    info['rootname'] = hdu[0].header['rootname']
    info['detector'] = segment
    info['temp'] = get_temp(full_path)

    times = timeline['time'][::step].copy()

    lat = timeline['latitude'][:-1][::step].copy().astype(np.float64)
    lon = timeline['longitude'][:-1][::step].copy().astype(np.float64)

    mjd = hdu[1].header['EXPSTART'] + \
        times.copy().astype(np.float64) * \
        SECOND_PER_MJD
    sun_lat = []
    sun_lon = []
    for item in get_sun_loc(mjd, full_path):
        sun_lon.append(item[0])
        sun_lat.append(item[1])

    mjd = mjd[:-1]

    decyear = mjd_to_decyear(mjd)

    if not len(times):
        logger.debug("time array empty for: {}".format(full_path))
        blank = np.array([0])
        yield info
        raise StopIteration

    events = hdu['events'].data
    filtered_index = np.where((events['PHA'] > pha[0]) &
                              (events['PHA'] < pha[1]) &
                              (events['XCORR'] > xlim[0]) &
                              (events['XCORR'] < xlim[1]) &
                              (events['YCORR'] > ylim[0]) &
                              (events['YCORR'] < ylim[1]))

    ta_index = np.where((events['XCORR'] > xlim[0]) &
                        (events['XCORR'] < xlim[1]) &
                        (events['YCORR'] > ylim[0]) &
                        (events['YCORR'] < ylim[1]))


    counts = np.histogram(events[filtered_index]['time'], bins=times)[0]
    ta_counts = np.histogram(events[ta_index]['time'], bins=times)[0]

    npix = float((xlim[1] - xlim[0]) * (ylim[1] - ylim[0]))
    counts = counts / npix / step
    ta_counts = ta_counts / npix / step

    if not len(lat) == len(counts):
        lat = lat[:-1]
        lon = lon[:-1]
        sun_lat = sun_lat[:-1]
        sun_lon = sun_lon[:-1]

    assert len(lat) == len(counts), \
        'Arrays are not equal in length {}:{}'.format(len(lat), len(counts))

    if not len(counts):
        logger.debug("zero-length array found for: {}".format(full_path))
        yield info
    else:
        for i in range(len(counts)):
            ### - better solution than round?
            info['date'] = round(decyear[i], 3)
            info['dark'] = round(counts[i], 7)
            info['ta_dark'] = round(ta_counts[i], 7)
            info['latitude'] = round(lat[i], 7)
            info['longitude'] = round(lon[i], 7)
            info['sun_lat'] = round(sun_lat[i], 7)
            info['sun_lon'] = round(sun_lon[i], 7)

            yield info

#-------------------------------------------------------------------------------

def compile_phd():
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

def make_plots(detector, base_dir, TA=False):
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
        solar_date = np.array( mjd_to_decyear([line[0] for line in solar_data]) )
        solar_flux = np.array([line[1] for line in solar_data])
    except TypeError:
        logger.warning("Couldn't read solar data.  Putting in all zeros.")
        solar_date = np.ones(1000)
        solar_flux = np.ones(1000)

    dark_key = 'dark'
    if TA:
        dark_key = 'ta_dark'

    #-- Open settings and get database
    settings = get_settings()
    database = get_database()

    for key, segment in zip(search_strings, segments):
        
        logger.debug('creating time plot for {}:{}'.format(segment, key))
        
        #-- Query for data here!
        data = Darks.select().where((Darks.detector == segment))
        
        data = [row for row in data]

        mjd = np.array([item.date for item in data])
        dark = np.array([item.dark for item in data])
        temp = np.array([item.temp for item in data])
        latitude = np.array([item.latitude for item in data])
        longitude = np.array([item.longitude for item in data])

        index = np.argsort(mjd)
        mjd = mjd[index]
        dark = dark[index]
        temp = temp[index]
        latitude = latitude[index]
        longitude = longitude[index]

        index_keep = np.where((longitude < 250) | (latitude > 10))[0]
        mjd = mjd[index_keep]
        dark = dark[index_keep]
        temp = temp[index_keep]

        outname = os.path.join(base_dir, detector, '{}_vs_time_{}.png'.format(dark_key, segment))
        if not os.path.exists(os.path.split(outname)[0]):
            os.makedirs(os.path.split(outname)[0])
        plot_time(detector, dark, mjd, temp, solar_flux, solar_date, outname)

        #-- Plot vs orbit
        logger.debug('creating orbit plot for {}:{}'.format(segment, key))
        
        data = Darks.select().where(Darks.detector==segment)
    
        data = [row for row in data]

        dark = np.array([item.dark for item in data])
        latitude = np.array([item.latitude for item in data])
        longitude = np.array([item.longitude for item in data])
        sun_lat = np.array([item.sun_lat for item in data])
        sun_lon = np.array([item.sun_lon for item in data])
        date = np.array([item.date for item in data])

        index = np.argsort(date)
        dark = dark[index]
        latitude = latitude[index]
        longitude = longitude[index]
        sun_lat = sun_lat[index]
        sun_lon = sun_lon[index]

        outname = os.path.join(base_dir, detector, '{}_vs_orbit_{}.png'.format(dark_key, segment))
        plot_orbital_rate(longitude, latitude, dark, sun_lon, sun_lat, outname)

        #-- Plot histogram of darkrates
        logger.debug('creating histogram plot for {}:{}'.format(segment, key))
        
        data = Darks.select().where(Darks.detector==segment)
        
        data = [item for item in data]

        dark = np.array([item.dark for item in data])
        date = np.array([item.date for item in data])

        index = np.argsort(date)
        date = date[index]
        dark = dark[index]

        for year in set(map(int, date)):

            if year == 0:
                continue
            else:
                index = np.where( (date >= year) &
                                (date < year + 1))
                outname = os.path.join(base_dir, detector, '{}_hist_{}_{}.png'.format(dark_key, year, segment))
                plot_histogram(dark[index], outname)
    
        index = np.where(date >= date.max() - .5)
        outname = os.path.join(base_dir, detector, '{}_hist_-6mo_{}.png'.format(dark_key, segment))
        plot_histogram(dark[index], outname )

        outname = os.path.join(base_dir, detector, '{}_hist_{}.png'.format(dark_key, segment))
        plot_histogram(dark, outname)

#-------------------------------------------------------------------------------

def move_products(base_dir, web_dir):
    '''Move created pdf files to webpage directory
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
                    logger.debug("removing {}".format(item))
                    move_list.remove(item)
                    continue
                else:
                    logger.debug("moving {}".format(item))
                
                #-- Split the file and paths.
                path, file_to_move = os.path.split(item)
                
                #-- Update the permissions
                os.chmod(item, 0o766)

                #-- Remove the file if it exists in the webpage dir.
                os.remove(write_dir + file_to_move)
                
                #-- Copy the file over.
                shutil.copy(item, write_dir + file_to_move)

            except OSError:
                logger.warning("Hit an os error for {}, leaving it there".format(item))
                move_list.remove(item)
        
        #-- Change all of the file permissions.
        os.system('chmod 777 ' + write_dir + '/*.png')
#-------------------------------------------------------------------------------

def monitor():
    """Main monitoring pipeline"""

    logger.info("Starting Monitor")

    settings = get_settings()
    out_dir = os.path.join(settings['monitor_location'], 'Darks')
    web_dir = settings['webpage_location']

    if not os.path.exists(out_dir):
        logger.warning("Creating output directory: {}".format(out_dir))
        os.makedirs(out_dir)

    get_solar_data(out_dir)

    for detector in ['FUV', 'NUV']:
        logger.info("Making plots for {}".format(detector))
        make_plots(detector, out_dir)

        if detector == 'FUV':
            make_plots(detector, out_dir, TA=True)

    logger.info("moving products to web directory")
    move_products(out_dir, web_dir)

#-------------------------------------------------------------------------------
