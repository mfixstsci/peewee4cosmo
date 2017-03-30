"""
Interactive plotting for Dark monitor
"""

from bokeh.io import output_file, show
from bokeh.layouts import column
from bokeh.plotting import figure

import scipy
from scipy.ndimage.filters import convolve
import numpy as np
import math
import os
from time import gmtime, strftime

#from ..utils import remove_if_there

#-------------------------------------------------------------------------------

def plot_time(detector, dark, date, temp, solar, solar_date, outname):
    """Plot the dar-rate vs time
    Parameters
    ----------
    detector : str
        FUV or NUV
    dark : np.ndarray
        array of measured dark rates in counts/s
    date : np.ndarray
        array of measured times
    temp : np.ndarray
        array of temperatures
    solar : np.ndarray
        array of solar flux values
    solar_date : np.ndarray
        array of solar dates
    outname : str
        path + name of output plot
    """
    
    #-- Check and then remove if file exists.
    #remove_if_there(outname)

    #-- Begin Python
    sorted_index = np.argsort(solar_date)
    solar = solar[sorted_index]
    solar_date = solar_date[sorted_index]

    #-- Plot height/width for FUV
    plt_hgt = 350
    plt_wth = 800

    #-- We track temperture for NUV mode
    if detector == 'NUV':
        temp_index = np.where(temp > 15)[0]
        dark = dark[temp_index]
        date = date[temp_index]
        temp = temp[temp_index]
        
        #-- More panels for NUV so make height shorter for browser...s
        plt_hgt = 250
        plt_wth = 800

    
    #-- Begin Bokeh   
    s1 = figure(width=plt_wth, height=plt_hgt, title='{} Global Dark Rate as of {}'.format(detector, strftime("%m-%d-%Y %H:%M:%S", gmtime())))
    s1.title.text_font_size = '15pt'
    s1.circle(date, dark, legend='Dark Count Rate',size=8, color="black", alpha=0.5)
    s1.yaxis.axis_label = "Mean Dark Rate (cnts/pix/sec)"

    #-- Plot NUV
    if detector  == 'NUV':

        #-- Smooth out solar data
        solar_smooth = scipy.convolve(solar, np.ones(81) / 81.0, mode='same')

        #-- Plot Solar Flux
        s2 = figure(width=plt_wth, height=plt_hgt, x_range=s1.x_range, title=None)
        s2.line(solar_date, solar, legend='10.7 cm', color='gold', line_width=2)
        s2.line(solar_date[:-41], solar_smooth[:-41], legend='10.7 cm smoothed', color='red', line_width=3)
        s2.yaxis.axis_label = "Radio Flux"
        s2.yaxis.axis_label_text_font_size = "15pt"

        #-- Plot Temperture
        s3 = figure(width=plt_wth, plot_height=plt_hgt, x_range=s1.x_range, title=None)
        s3.circle(date, temp, size=8, color="red", alpha=0.5)
        s3.xaxis.axis_label = "Decimal Year"
        s3.yaxis.axis_label = "Temperture"
        
        s3.xaxis.axis_label_text_font_size = "15pt"
        s3.yaxis.axis_label_text_font_size = "15pt"

        p = column(s1, s2, s3)

    #-- Plot FUV
    if detector  == 'FUV':

        #-- FUV only has two panels, make the font bigger
        s1.yaxis.axis_label_text_font_size = "15pt"

        #-- Smooth out solar data
        solar_smooth = scipy.convolve(solar, np.ones(81) / 81.0, mode='same')

        #-- Plot Solar Flux
        s2 = figure(width=plt_wth, height=plt_hgt, x_range=s1.x_range, title=None)
        s2.line(solar_date, solar, legend='10.7 cm', color='gold', line_width=2)
        s2.line(solar_date[:-41], solar_smooth[:-41], legend='10.7 cm smoothed', color='red', line_width=3)
        s2.yaxis.axis_label = "Radio Flux"
        s2.xaxis.axis_label = "Decimal Year"
        s2.xaxis.axis_label_text_font_size = "15pt"
        s2.yaxis.axis_label_text_font_size = "15pt"

        p = column(s1, s2)

    #-- Save file
    output_file(outname)
#-------------------------------------------------------------------------------
