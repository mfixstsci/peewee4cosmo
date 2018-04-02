"""Create periodogram of COS FUV dark rate for a given time period.
Trying to constrain the cause in the variation of the dark rate
based off of periodogram.

"""

__author__ = 'Mees Fix'
__maintainer__ = 'Mees Fix'
__email__ = 'mfix@stsci.edu'
__status__ = 'Active'

from astropy.stats import LombScargle
from astropy.time import Time

import matplotlib.pyplot as plt  

from peewee import *

from ..database.models import get_database, get_settings
from ..database.models import Darks

import os

from bokeh.io import output_file, show, save
from bokeh.plotting import figure

import numpy as np

from scipy.fftpack import fft


def create_periodogram(x, y, outname):
    """Create periodogram based off of astropy's LombScagle method in stats
    library.

    Parameters
    ----------
    x: list like
        Variable in x (dark = Time)
    y: list like 
        Variable in y (dark = cnts/pix/sec)
    segment: str
        FUVA or FUVB
    Creates
    -------
    periodogram
    """
    
    plt_hgt = 600
    plt_wth = 900

    output_file(outname)

    frequency, power = LombScargle(x, y).autopower()
    
    p = figure(width=plt_wth, height=plt_hgt, title='Periodogram')
    p.title.text_font_size = '15pt'
    p.line(frequency, power, color="black")
    p.xaxis.axis_label = "Frequency"
    p.yaxis.axis_label = "Power"

    save(p, filename=outname)


def create_periodogram(x, y, outname):
    """Create periodogram based off of astropy's LombScagle method in stats
    library.

    Parameters
    ----------
    x: list like
        Variable in x (dark = Time)
    y: list like 
        Variable in y (dark = cnts/pix/sec)
    segment: str
        FUVA or FUVB
    Creates
    -------
    periodogram
    """
    
    plt_hgt = 600
    plt_wth = 900

    output_file(outname)

    frequency, power = LombScargle(x, y).autopower()
    
    p = figure(width=plt_wth, height=plt_hgt, title='Periodogram')
    p.title.text_font_size = '15pt'
    p.line(frequency, power, color="black")
    p.xaxis.axis_label = "Frequency"
    p.yaxis.axis_label = "Power"

    save(p, filename=outname)


def dark_query(start_date, end_date, segment):
    """Query the dark data in COSMO between two dates.

    Parameters
    ----------
    start_date: float
        Beginning date of dark data
    end_date: float
        Date end point for dark data.
    segment: str
        FUVA or FUVB
    """

    if start_date > end_date:
        start_date, end_date = end_date, start_date
        print('DATES ARE FLIPPED!')
    
    database = get_database()
    database.connect()

    data = Darks.select().order_by(Darks.date)\
                            .where(
                                   (Darks.date >= start_date)
                                   & (Darks.date <= end_date)
                                   & (Darks.detector == segment))
     
    date = Time([row.date for row in data], format='jyear').seconds
    dark = [row.dark for row in data]
    print(date.max()-date.min())
    return date, dark 


def main():
    """Main driver.
    """
    settings = get_settings()
    for segment in ['FUVA', 'FUVB']:
        out_dir = os.path.join(settings['monitor_location'], 'Darks', 'FUV', 
                               'dark_periodicty_{}.html'.format(segment))
        date, dark = dark_query(2016, 2018, segment)
        create_periodogram(date, dark, out_dir)