"""
Trend gain for pixel-to-pixel changes and store pixel as row in DB Table. 
Row for each pixel is updated each time scripts are run.
"""

import os
import numpy as np

from astropy.io import fits

import scipy
from scipy.optimize import leastsq, newton, curve_fit

from ..database.models import get_database, get_settings
from ..database.models import Gain_Trends, Gain
from .constants import * 

import collections
import functools
import multiprocessing as mp

import scipy
from scipy.optimize import leastsq, newton, curve_fit

import itertools

import logging
logger = logging.getLogger(__name__)
#-------------------------------------------------------------------------------

def time_trends(segment, hv_lvl, out_dir):

    database = get_database()
    database.connect()

    results = Gain_Trends.select().where(
                                         (Gain_Trends.segment == segment) &
                                         (Gain_Trends.hv_lvl == hv_lvl)  
                                         )

    slope_image = np.zeros((1024, 16384))
    intercept_image = np.zeros((1024, 16384))
    bad_image = np.zeros((1024, 16384))

    for row in results:
        y = row.y * Y_BINNING
        x = row.x * X_BINNING
        slope_image[y:y+Y_BINNING, x:x+X_BINNING] = row.slope
        intercept_image[y:y+Y_BINNING, x:x+X_BINNING] = row.intercept

        bad_image[y:y+Y_BINNING, x:x+X_BINNING] = row.mjd

    if slope_image.any():
        logger.debug("WRITING PROJECTION FOR {} {}".format(segment, hv_lvl))
        write_projection(out_dir, slope_image, intercept_image, bad_image, segment, hv_level)

#-------------------------------------------------------------------------------
def measure_slopes(coords, gain_dict, segment, hv_lvl):
    """
    Calculate the slope for pixel
    """

    data = {}
    x, y = coords
    all_gain = np.array([gain['gain'] for gain in gain_dict[coords]])
    all_expstart = np.array([expstart['expstart'] for expstart in gain_dict[coords]])

    data['segment'] = segment
    data['hv_lvl'] = hv_lvl
    data['x'] = x
    data['y'] = y

    #-- Want at least 5 measurements for the fit.
    if not len(all_gain) > 5:
        data['slope'] = 0.0
        data['intercept'] = 0.0
        data['proj_bad_mjd'] = 0.0
        return data

    #-- Sort indices based on gain
    sorted_index = all_gain.argsort()
    all_gain = all_gain[sorted_index]
    all_expstart = all_expstart[sorted_index]

    #-- Fit line to make prediciton
    fit, parameters, success = time_fitting(all_expstart, all_gain)

    #-- If the fit was successful, add to db.
    if success:
        intercept = parameters[1]
        slope = parameters[0]

        f = lambda x, a, b: a * x + b - 3
        fprime = lambda x, a, b: a

        try:
            date_bad = newton(f, all_gain[-1], fprime, args=tuple(parameters), tol=1e-5, maxiter=1000)
        except RuntimeError:
            date_bad = 0.0
        
        
        data['slope'] = round(slope, 5)
        data['intercept'] = round(intercept, 5)
        data['proj_bad_mjd'] = round(date_bad, 5)

    return data
#-------------------------------------------------------------------------------
def query_and_sort_gain(segment, hv_level):
    """
    Measure slope for each pixel.
    """

    #-- Connect to database.
    database = get_database()
    database.connect()
    settings = get_settings()

    flagged_pix = list(Gain.select().distinct().where(
                                                      (Gain.segment == segment) &
                                                      (Gain.hv_lvl == hv_level) &
                                                      (Gain.y.between(400//Y_BINNING, 600//Y_BINNING))
                                                      )
                                                     .dicts())
    
    #-- Create list of dictionaries where the key is the x,y pair.
    logger.info('ADDING NEW TREND ENTRIES FOR SEGMENT {}, HV {}'.format(segment, hv_level))
    result = collections.defaultdict(list)
    for d in flagged_pix:
        result[d['x'], d['y']].append(d)
    
    partial = functools.partial(measure_slopes,
                                segment = segment,
                                gain_dict=result,
                                hv_lvl=hv_level)
    
    pool = mp.Pool(processes=settings['num_cpu'])                        
    trend_data = pool.map(partial, result.keys())

    with database.atomic():
        Gain_Trends.insert_many(trend_data).execute() 
    
    database.close()

#-------------------------------------------------------------------------------

def time_fitting(x_fit, y_fit):
    """Fit a linear relation to the x_fit and y_fit parameters
    Parameters
    ----------
    x_fit : np.ndarray
        x-values to fit
    y_fit : np.ndarray
        y-values to fit
    Returns
    -------
    fit, parameters, success : tuple
        fit to the values, fit parameters, boolean-success
    """

    x_fit = np.array(x_fit)
    y_fit = np.array(y_fit)

    ###First fit iteration and remove outliers
    POLY_FIT_ORDER = 1

    slope, intercept = scipy.polyfit(x_fit, y_fit, POLY_FIT_ORDER)
    fit = scipy.polyval((slope, intercept), x_fit)
    fit_sigma = fit.std()
    include_index = np.where(np.abs(fit-y_fit) < 1.5*fit_sigma)[0]

    if len(include_index) < 4:
        return None, None, False

    x_fit_clipped = x_fit[include_index]
    y_fit_clipped = y_fit[include_index]

    parameters = scipy.polyfit(x_fit_clipped, y_fit_clipped, POLY_FIT_ORDER)
    fit = scipy.polyval(parameters, x_fit)

    return fit, parameters, True

#-------------------------------------------------------------------------------

def write_projection(out_dir, slope_image, intercept_image, bad_image, segment, dethv):
    """Writs a fits file with information useful for post-monitoring analysis.
    Parameters
    ----------
    slope_image : np.ndarray
        2D image of linear gain degredation slopes
    intercept_image : np.ndarray
        2D image of intercepts for the linear gain degredations
    bad_image : np.ndarray
        2D image of the extrapolated date where the gain will drop below 3
    segment : str
        'FUVA' or 'FUVB', COS detector segment of the measurements
    dethv : int
        Detector high-voltage setting of the measurements
    Returns
    -------
        None
    Outputs
    -------
        FITS file with the saved array data.
    """

    hdu_out = fits.HDUList(fits.PrimaryHDU())
    hdu_out[0].header['TELESCOP'] = 'HST'
    hdu_out[0].header['INSTRUME'] = 'COS'
    hdu_out[0].header['DETECTOR'] = 'FUV'
    hdu_out[0].header['OPT_ELEM'] = 'ANY'
    hdu_out[0].header['Fsave_dirILETYPE'] = 'PROJ_BAD'
    hdu_out[0].header['DETHV'] = dethv
    hdu_out[0].header['SEGMENT'] = segment

    #---Ext 1
    hdu_out.append(fits.ImageHDU(data=bad_image))
    hdu_out[1].header['EXTNAME'] = 'PROJBAD'

    #---Ext 2
    hdu_out.append(fits.ImageHDU(data=slope_image))
    hdu_out[2].header['EXTNAME'] = 'SLOPE'

    #---Ext 3
    hdu_out.append(fits.ImageHDU(data=intercept_image))
    hdu_out[3].header['EXTNAME'] = 'INTERCEPT'

    #---Writeout
    hdu_out.writeto(os.path.join(out_dir, 'proj_bad_{}_{}.fits'.format(segment, dethv)) ,clobber=True)
    hdu_out.close()

#-------------------------------------------------------------------------------

def main():
    
    settings = get_settings()
    database = get_database()
    database.connect()

    all_combos = [(combo.segment, combo.hv_lvl) for combo in \
                  Gain.select(Gain.segment, Gain.hv_lvl).distinct()]
    database.close()

    for segment, hv_level in all_combos:
        query_and_sort_gain(segment, hv_level)
    
    out_dir = os.path.join(settings['monitor_location'], 'CCI')
    
    for segment, hv_level in all_combos:
        time_trends(segment, hv_level, out_dir)
    