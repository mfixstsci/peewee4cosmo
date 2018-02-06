from __future__ import absolute_import

""" Script to compile the spectrum shift data for COS FUV and NUV data.
"""

import glob
import os
import shutil
import sys
import logging
logger = logging.getLogger(__name__)

import numpy as np
np.seterr(divide='ignore')
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io import ascii

import scipy
from scipy.stats import linregress
from datetime import datetime

from astropy.io import fits
from astropy.table import Table
from astropy import table
from time import gmtime, strftime, localtime

from ..database.models import get_database, get_settings
from ..database.models import Lampflash, Rawacqs, Files

from ..utils import remove_if_there

from copy import deepcopy

import pandas as pd

from bokeh.io import output_file, show, save
from bokeh.plotting import figure
from bokeh.models import Range1d, HoverTool, BoxSelectTool, ColumnDataSource, OpenURL, TapTool, Div, Button, CustomJS
from bokeh.layouts import column, row
#-------------------------------------------------------------------------------

def fppos_shift(lamptab_name, segment, opt_elem, cenwave, fpoffset):
    """Get the COS FPPOS pixel shift.
    
    Parameters
    ----------
    lamptab_name : str
        name of lamptab reference file.
    segment : str
        detector segment.
    opt_elem : str
        grating
    cenwave : str
        central wavelength
    fpoffset : int
        calculated offset from home position.
    """

    lamptab = fits.getdata(os.path.join(os.environ['lref'], lamptab_name))

    if 'FPOFFSET' not in lamptab.names:
        return 0

    index = np.where((lamptab['segment'] == segment) &
                     (lamptab['opt_elem'] == opt_elem) &
                     (lamptab['cenwave'] == cenwave) &
                     (lamptab['fpoffset'] == fpoffset))[0]

    offset = lamptab['FP_PIXEL_SHIFT'][index][0]

    return offset

#-------------------------------------------------------------------------------

def pull_flashes(filename):
    """Calculate lampflash values for given file

    Parameters
    ----------
    filename : str
        file to calculate lamp shifts from

    Returns
    -------
    out_info : dict
        dictionary of pertinent value

    """

    #-- Open file
    file_path = os.path.join(filename.path, filename.filename)
    with fits.open(file_path) as hdu:
        #-- Set some dictionary values.
        out_info = {'filename': filename.filename,
                    'date': hdu[1].header['EXPSTART'],
                    'rootname': hdu[0].header['ROOTNAME'],
                    'proposid': hdu[0].header['PROPOSID'],
                    'detector': hdu[0].header['DETECTOR'],
                    'opt_elem': hdu[0].header['OPT_ELEM'],
                    'cenwave': hdu[0].header['CENWAVE'],
                    'fppos': hdu[0].header.get('FPPOS', None),
                    'filetype': hdu[0].header.get('FILETYPE', None)}

        #-- Get time, and then convert the format
        t = Time(out_info['date'], format='mjd')
        out_info['cal_date'] = t.iso

        #-- Open lampflash
        if '_lampflash.fits' in filename.filename:
            
            out_info['segment'] = hdu[0].header['SEGMENT']
            #-- Get lamptab file
            out_info['lamptab'] = hdu[0].header['LAMPTAB'].split('$')[-1]

            #-- FPPOS 3 is the home frame, so put all FP's in home frame.
            fpoffset = out_info['fppos'] - 3

            for i, line in enumerate(hdu[1].data):
                #-- Count the number of flashes and set dictionary values.
                out_info['flash'] = (i // 2) + 1
                out_info['x_shift'] = line['SHIFT_DISP'] - fppos_shift(out_info['lamptab'],
                                                                        line['segment'],
                                                                        out_info['opt_elem'],
                                                                        out_info['cenwave'],
                                                                        fpoffset)

                out_info['y_shift'] = line['SHIFT_XDISP']
                out_info['found'] = line['SPEC_FOUND']
                out_info['segment'] = line['SEGMENT']

                #-- don't need too much precision here
                out_info['x_shift'] = round(out_info['x_shift'], 5)
                out_info['y_shift'] = round(out_info['y_shift'], 5)

                yield deepcopy(out_info)
        
        #-- Open rawacqs
        elif '_rawacq.fits' in filename.filename:
            #-- Technically it wasn't found.
            out_info['found'] = False
            out_info['fppos'] = -1
            out_info['flash'] = 1

            #-- Grab associated spt
            spt = fits.open(os.path.join(filename.path,filename.filename.replace('rawacq', 'spt')))
            
            if not spt[1].header['LQTAYCOR'] > 0:
                out_info['x_shift'] = -999
                out_info['y_shift'] = -999
            else:
                # These are in COS RAW coordinates, so shifted 90 degrees from
                # user and backwards
                out_info['x_shift'] = 1023 - spt[1].header['LQTAYCOR']
                out_info['y_shift'] = 1023 - spt[1].header['LQTAXCOR']

            yield deepcopy(out_info)
        else:
            yield deepcopy(out_info)

#-------------------------------------------------------------------------------

def fit_data(xdata, ydata):
    """ Fit a regression line to shift data points

    Parameters
    ----------
    xdata : astropy.table.column.Column
        A list of x values (time)
    ydata : astropy.table.column.Column
        A list of y values (shifts)

    Returns
    -------
    fit : ndarray
        The fit line
    xdata : astropy.table.column.Column
        List of x values for fit
    parameters : tuple
        fitting parameters
    err : int
        Value returned on whether the fit was a sucess.
    """
    stats = linregress(xdata, ydata)

    parameters = (stats[0], stats[1])
    err = 0
    fit = scipy.polyval(parameters, xdata)

    return fit, xdata, parameters, err

#-------------------------------------------------------------------------------

def make_shift_table(db_table):
    """ Make an astropy table of shift values and other metadata
    
    Parameters
    ----------
    db_table : peewee table object
        The Lampflash or Rawacq table
    
    Returns
    -------
    data : Astropy table
        All data needed for plotting obtained from database.
    
    """
    database = get_database()
    database.connect()

    data = []

    #-- this is a crude implementation, but it lets me use the rest of the
    #-- plotting code as-is

    #-- .dicts() returns the result objects as dictionaries. 
    for i, row in enumerate(db_table.select().dicts()):
        data.append(row.values())
        if not i:
            #-- get keys here because if you use ._meta.fields.keys() 
            #-- they will be out of order.
            keys = row.keys()    
    
    database.close()

    data = Table(rows=data, names=keys)

    return data

#-------------------------------------------------------------------------------

def make_panel(data, grating, height, width, detector, plt_color, top=False, x_range=False, acqs=False, all_nuv=False):
    """Make a bokeh panel for figure.

    Parameters
    ----------
    data: Astropy.Table
        Astropy table of all metadata.
    grating: np.array
        Data associated with a COS grating from table.
    height: int
        Height of panel.
    width: int
        Width of panel.
    detector: str
        FUV or NUV.
    plt_color: str
        color of scatter points.
    top: Bool
        Is it the top panel of the figure, True or False.
    x_range: 
    """

    #-- Define tools that each panel will possess.
    TOOLS ='box_zoom,box_select,pan,reset,tap'
    
    #-- Build ColumnDataSource object
    source = ColumnDataSource(data=dict(
                    date=data['date'][grating],
                    shift=data['x_shift'][grating],
                    proposid=data['proposid'][grating],
                    rootname=data['rootname'][grating],
                    ))

    hover = HoverTool(
                tooltips=[
                        ("Time", "@date"),
                        ("Shift", "@shift"),
                        ("Proposid", "@proposid"),
                        ("Rootname", "@rootname"),
                        ]
                    )

    #-- Parse detector for labeling
    if detector == 'FUV':
        #-- If top panel, add a title
        if top:
            panel = figure(width=width, height=height, x_range=(min(data['date']) - 10, max(data['date']) + 10), title='FUV SHIFT1[A/B] as of {} EST'.format(strftime("%m-%d-%Y %H:%M:%S", localtime())), tools=[TOOLS,hover])
            panel.title.text_font_size = '15pt'
        else:
            panel = figure(width=width, height=height, x_range=x_range, title=None, tools=[TOOLS, hover])
       
        #-- Label y and also draw max bounds from reference file. 
        panel.line(data['date'], np.zeros_like(data['date']) + 300, color='black', line_width=2, line_dash='dashed')
        panel.line(data['date'], np.zeros_like(data['date']) - 300, color='black', line_width=2, line_dash='dashed')
        panel.yaxis.axis_label = "Shift1[A/B] (Pixels)"
    
    elif detector == 'NUV':
        #-- If top panel, add a title
        if top:
            panel = figure(width=width, height=height, x_range=(min(data['date']) - 10, max(data['date']) + 10), title='NUV SHIFT1[A/B/C] as of {} EST'.format(strftime("%m-%d-%Y %H:%M:%S", localtime())), tools=[TOOLS, hover])
            panel.title.text_font_size = '15pt'            
        else:
            panel = figure(width=width, height=height, x_range=x_range, title=None, tools=[TOOLS, hover])
        panel.yaxis.axis_label = "Shift1[A/B/C] (Pixels)"

    #-- Make scatter plot of data
    if all_nuv:
        panel.circle('date', 'shift', legend='All NUV', size=4, source=source, color=plt_color, alpha=0.5)
    else:
        panel.circle('date', 'shift', legend=data['opt_elem'][grating][0], size=4, source=source, color=plt_color, alpha=0.5)

    #-- Provide URL and taptool and callback info.
    url = "http://archive.stsci.edu/proposal_search.php?id=@proposid&mission=hst"
    taptool = panel.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    
    if not acqs:
        #-- Draw a line @ 0.
        panel.line(data['date'][grating], np.zeros_like(data['date'][grating]), color='red', line_width=2)
    
    return panel
#-------------------------------------------------------------------------------

def make_interactive_plots(data, data_acqs, out_dir, detector):
    """Make interactive plots for OSM shifts  
    
    Parameter
    ---------
    data : Astropy Table
        A table of lampflash metadata
    data_acqs : Astropy Table
        A table of rawacqs metadata
    out_dir : str
        The output directory for the files.
    detector : str
        FUV or NUV mode to make correct plot.
    """

    logger.info("MAKING INTERACTIVE PLOT FOR {}".format(detector))

    #-- Sort by time
    sorted_index = np.argsort(data['date'])
    
    #-- Rearange Data
    data = data[sorted_index]

    if detector == 'FUV':
        
        #-- Sort all data by opt_elem
        G140L = np.where((data['opt_elem'] == 'G140L'))[0]
        G130M = np.where((data['opt_elem'] == 'G130M'))[0]
        G160M = np.where((data['opt_elem'] == 'G160M'))[0]

        #-- Find Unique Entries
        unique_data = table.unique(data, ['date','x_shift','rootname','proposid'])
        unique_G140L = np.where((unique_data['opt_elem'] == 'G140L'))[0]
        unique_G130M = np.where((unique_data['opt_elem'] == 'G130M'))[0]
        unique_G160M = np.where((unique_data['opt_elem'] == 'G160M'))[0]

        #-- Begin Bokeh   
        TOOLS ='box_zoom,box_select,crosshair,pan,reset,tap'
        
        #-- Plot FUV Shifts
        outname = os.path.join(out_dir, 'FUV_shift_vs_time.html')
        remove_if_there(outname)
        output_file(outname)
        
        #-- Set panel size
        plt_hgt = 250
        plt_wth = 800

        #-- Create bokeh figure objects
        
        #-- Panel 1
        s1 = make_panel(unique_data, unique_G130M, plt_hgt, plt_wth, 'FUV', 'blue', top=True)
        #-- Fit shift as a function of date and plot it...
        fit,ydata,parameters,err = fit_data(data['date'][G130M],data['x_shift'][G130M])
        s1.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))
        
        #-- Panel 2
        s2 = make_panel(unique_data, unique_G160M, plt_hgt, plt_wth, 'FUV', 'green', x_range=s1.x_range)
        #-- Fit shift as a function of date and plot it...
        fit,ydata,parameters,err = fit_data(data['date'][G160M],data['x_shift'][G160M])
        s2.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))
        
        #-- Panel 3
        s3 = make_panel(unique_data, unique_G140L, plt_hgt, plt_wth, 'FUV', 'yellow', x_range=s1.x_range)
        
        #-- Fit shift as a function of date and plot it...
        fit,ydata,parameters,err = fit_data(data['date'][G140L],data['x_shift'][G140L])
        s3.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))
        
        s3.xaxis.axis_label = "Time (MJD)"
        
        p = column(s1, s2, s3)

        save(p, filename=outname)

        return p
        
    #-- NUV Plots
    if detector == 'NUV':
        
        #-- Sort by grating.
        G230L = np.where((data['opt_elem'] == 'G230L'))[0]
        G225M = np.where((data['opt_elem'] == 'G225M'))[0]
        G285M = np.where((data['opt_elem'] == 'G285M'))[0]
        G185M = np.where((data['opt_elem'] == 'G185M'))[0]
        NUV = np.where((data['opt_elem'] == 'G230L') |
                    (data['opt_elem'] == 'G185M') |
                    (data['opt_elem'] == 'G225M') |
                    (data['opt_elem'] == 'G285M'))[0]
        mirrora = np.where((data_acqs['opt_elem'] == 'MIRRORA')
                        & (data_acqs['x_shift'] > 0))[0]
        mirrorb = np.where((data_acqs['opt_elem'] == 'MIRRORB')
                        & (data_acqs['x_shift'] > 0))[0]


        #-- Sort unique entries by grating.
        unique_data = table.unique(data, ['date','x_shift','rootname','proposid'])
        unique_acqs = table.unique(data_acqs, ['date','x_shift','rootname','proposid'])

        unique_G230L = np.where((unique_data['opt_elem'] == 'G230L'))[0]
        unique_G225M = np.where((unique_data['opt_elem'] == 'G225M'))[0]
        unique_G285M = np.where((unique_data['opt_elem'] == 'G285M'))[0]
        unique_G185M = np.where((unique_data['opt_elem'] == 'G185M'))[0]
        unique_NUV = np.where((unique_data['opt_elem'] == 'G230L') |
                    (unique_data['opt_elem'] == 'G185M') |
                    (unique_data['opt_elem'] == 'G225M') |
                    (unique_data['opt_elem'] == 'G285M'))[0]
        unique_mirrora = np.where((unique_acqs['opt_elem'] == 'MIRRORA')
                        & (unique_acqs['x_shift'] > 0))[0]
        unique_mirrorb = np.where((unique_acqs['opt_elem'] == 'MIRRORB')
                        & (unique_acqs['x_shift'] > 0))[0]
        
        #-- Bokeh
        TOOLS ='box_zoom,pan,reset,hover,tap'

        #-- Bokeh panel sizes.
        plt_hgt = 250
        plt_wth = 800

        #-- Set outname and create file.
        outname = os.path.join(out_dir, 'NUV_shift_vs_time.html')
        remove_if_there(outname)
        output_file(outname)
        
        #-- G230L search range was updated earlier than the other observing modes.
        transition_date = 56500.0
        transition_date_G230L = 55535.0        

        #-- Because of the complexity of the different transition dates, plotting with bokeh and acq figures... I've made code blocks for each panel. 
        #-- Panel 1
        #-- Create bokeh figure.
        s1 = make_panel(unique_data, unique_G185M, plt_hgt, plt_wth, 'NUV', 'blue', top=True)
        #-- Fit Data
        fit,ydata,parameters,err = fit_data(data['date'][G185M],data['x_shift'][G185M])
        s1.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))

        #-- Find transition regions.
        before_data = np.where(data['date'][G185M] <= transition_date)
        after_data = np.where(data['date'][G185M] >= transition_date)
        
        #-- First transitions
        s1.line(data['date'][G185M][before_data], np.zeros_like(data['date'][G185M][before_data]) + 58, color='black', line_width=2, line_dash='dashed')
        s1.line(data['date'][G185M][before_data], np.zeros_like(data['date'][G185M][before_data]) - 58, color='black', line_width=2, line_dash='dashed')

        #-- Second
        s1.line(data['date'][G185M][after_data], np.zeros_like(data['date'][G185M][after_data]) + 90, color='black', line_width=2, line_dash='dashed')
        s1.line(data['date'][G185M][after_data], np.zeros_like(data['date'][G185M][after_data]) - 90, color='black', line_width=2, line_dash='dashed')
        ############################

        #-- Panel 2
        s2 = make_panel(unique_data, unique_G225M, plt_hgt, plt_wth, 'NUV', 'red', x_range=s1.x_range)

        fit,ydata,parameters,err = fit_data(data['date'][G225M],data['x_shift'][G225M])
        s2.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))
        before_data = np.where(data['date'][G225M] <= transition_date)
        after_data = np.where(data['date'][G225M] >= transition_date)
        
        s2.line(data['date'][G225M][before_data], np.zeros_like(data['date'][G225M][before_data]) + 58, color='black', line_width=2, line_dash='dashed')
        s2.line(data['date'][G225M][before_data], np.zeros_like(data['date'][G225M][before_data]) - 58, color='black', line_width=2, line_dash='dashed')

        s2.line(data['date'][G225M][after_data], np.zeros_like(data['date'][G225M][after_data]) + 90, color='black', line_width=2, line_dash='dashed')
        s2.line(data['date'][G225M][after_data], np.zeros_like(data['date'][G225M][after_data]) - 90, color='black', line_width=2, line_dash='dashed')
        ############################

        #-- Panel 3
        s3 = make_panel(unique_data, unique_G285M, plt_hgt, plt_wth, 'NUV', 'yellow', x_range=s1.x_range)
        
        fit,ydata,parameters,err = fit_data(data['date'][G285M],data['x_shift'][G285M])
        s3.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))
        
        before_data = np.where(data['date'][G285M] <= transition_date)
        after_data = np.where(data['date'][G285M] >= transition_date)
        s3.line(data['date'][G285M][before_data], np.zeros_like(data['date'][G285M][before_data]) + 58, color='black', line_width=2, line_dash='dashed')
        s3.line(data['date'][G285M][before_data], np.zeros_like(data['date'][G285M][before_data]) - 58, color='black', line_width=2, line_dash='dashed')

        s3.line(data['date'][G285M][after_data], np.zeros_like(data['date'][G285M][after_data]) + 90, color='black', line_width=2, line_dash='dashed')
        s3.line(data['date'][G285M][after_data], np.zeros_like(data['date'][G285M][after_data]) - 90, color='black', line_width=2, line_dash='dashed')
        ############################
        
        #-- Panel 4
        s4 = make_panel(unique_data, unique_G230L, plt_hgt, plt_wth, 'NUV', 'green', x_range=s1.x_range)

        fit,ydata,parameters,err = fit_data(data['date'][G230L],data['x_shift'][G230L])
        s4.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))

        before_data = np.where(data['date'][G230L] <= transition_date_G230L)
        after_data = np.where(data['date'][G230L] >= transition_date_G230L)
        s4.line(data['date'][G230L][before_data], np.zeros_like(data['date'][G230L][before_data]) + 58, color='black', line_width=2, line_dash='dashed')
        s4.line(data['date'][G230L][before_data], np.zeros_like(data['date'][G230L][before_data]) - 58, color='black', line_width=2, line_dash='dashed')

        s4.line(data['date'][G230L][after_data], np.zeros_like(data['date'][G230L][after_data]) + 90, color='black', line_width=2, line_dash='dashed')
        s4.line(data['date'][G230L][after_data], np.zeros_like(data['date'][G230L][after_data]) - 90, color='black', line_width=2, line_dash='dashed')
        ############################
        
        #-- Panel 5
        s5 = make_panel(unique_data, unique_NUV, plt_hgt, plt_wth, 'NUV', 'firebrick', x_range=s1.x_range, all_nuv=True)
        fit,ydata,parameters,err = fit_data(data['date'][NUV],data['x_shift'][NUV])
        s5.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))
        ############################
        
        #-- Panel 6
        s6 = make_panel(unique_acqs, unique_mirrora, plt_hgt, plt_wth, 'NUV', 'firebrick', x_range=s1.x_range, acqs=True)
        fit,ydata,parameters,err = fit_data(data_acqs['date'][mirrora],data_acqs['x_shift'][mirrora])
        s6.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))
        ############################
        

        #-- Panel 7
        s7 = make_panel(unique_acqs, unique_mirrorb, plt_hgt, plt_wth, 'NUV', 'firebrick', x_range=s1.x_range, acqs=True)
        fit,ydata,parameters,err = fit_data(data_acqs['date'][mirrorb],data_acqs['x_shift'][mirrorb])
        s7.line(ydata, fit, color='black', line_width=2, legend=str(parameters[0]))
        s7.xaxis.axis_label = "Date (MJD)"
        ############################

        #-- Format into single column.
        p = column(s1, s2, s3, s4, s5, s6, s7)

        save(p, filename=outname)

        return p        
#-------------------------------------------------------------------------------

def make_plots(data, data_acqs, out_dir):
    """Make plots for OSM shifts  
    
    Parameter
    ---------
    data : Astropy Table
        A table of lampflash metadata
    data_acqs : Astropy Table
        A table of rawacqs metadata
    out_dir : str
        The output directory for the files.
    """
    
    logger.info("MAKING STATIC PLOTS")

    mpl.rcParams['figure.subplot.hspace'] = 0.05
    
    plt.rc('font', weight='bold')
    plt.rc('xtick.major', size=5, pad=7)
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=10)   
    
    sorted_index = np.argsort(data['date'])
    data = data[sorted_index]

    G140L = np.where((data['opt_elem'] == 'G140L'))[0]
    G140L_A = np.where((data['opt_elem'] == 'G140L') &
                       (data['segment'] == 'FUVA'))[0]
    G140L_B = np.where((data['opt_elem'] == 'G140L') &
                       (data['segment'] == 'FUVB'))[0]

    G130M = np.where((data['opt_elem'] == 'G130M'))[0]
    G130M_A = np.where((data['opt_elem'] == 'G130M') &
                       (data['segment'] == 'FUVA'))[0]
    G130M_B = np.where((data['opt_elem'] == 'G130M') &
                       (data['segment'] == 'FUVB'))[0]

    G160M = np.where((data['opt_elem'] == 'G160M'))[0]
    G160M_A = np.where((data['opt_elem'] == 'G160M') &
                       (data['segment'] == 'FUVA'))[0]
    G160M_B = np.where((data['opt_elem'] == 'G160M') &
                       (data['segment'] == 'FUVB'))[0]

    G230L = np.where((data['opt_elem'] == 'G230L'))[0]
    G230L_A = np.where((data['opt_elem'] == 'G230L') &
                       (data['segment'] == 'NUVA'))[0]
    G230L_B = np.where((data['opt_elem'] == 'G230L') &
                       (data['segment'] == 'NUVB'))[0]
    G230L_C = np.where((data['opt_elem'] == 'G230L') &
                       (data['segment'] == 'NUVC'))[0]

    G225M = np.where((data['opt_elem'] == 'G225M'))[0]
    G225M_A = np.where((data['opt_elem'] == 'G225M') &
                       (data['segment'] == 'NUVA'))[0]
    G225M_B = np.where((data['opt_elem'] == 'G225M') &
                       (data['segment'] == 'NUVB'))[0]
    G225M_C = np.where((data['opt_elem'] == 'G225M') &
                       (data['segment'] == 'NUVC'))[0]

    G285M = np.where((data['opt_elem'] == 'G285M'))[0]
    G285M_A = np.where((data['opt_elem'] == 'G285M') &
                       (data['segment'] == 'NUVA'))[0]
    G285M_B = np.where((data['opt_elem'] == 'G285M') &
                       (data['segment'] == 'NUVB'))[0]
    G285M_C = np.where((data['opt_elem'] == 'G285M') &
                       (data['segment'] == 'NUVC'))[0]

    G185M = np.where((data['opt_elem'] == 'G185M'))[0]
    G185M_A = np.where((data['opt_elem'] == 'G185M') &
                       (data['segment'] == 'NUVA'))[0]
    G185M_B = np.where((data['opt_elem'] == 'G185M') &
                       (data['segment'] == 'NUVB'))[0]
    G185M_C = np.where((data['opt_elem'] == 'G185M') &
                       (data['segment'] == 'NUVC'))[0]

    NUV = np.where((data['opt_elem'] == 'G230L') |
                   (data['opt_elem'] == 'G185M') |
                   (data['opt_elem'] == 'G225M') |
                   (data['opt_elem'] == 'G285M'))[0]

    #############

    fig = plt.figure( figsize=(16,8) )
    ax = fig.add_subplot(3,1,1)

    ax.plot( data['date'][G130M_A], data['x_shift'][G130M_A],'b.',label='G130M')
    ax.plot( data['date'][G130M_B], data['x_shift'][G130M_B],'b.')
    ax.xaxis.set_ticklabels( ['' for item in ax.xaxis.get_ticklabels()] )

    ax2 = fig.add_subplot(3,1,2)
    ax2.plot( data['date'][G160M_A], data['x_shift'][G160M_A],'g.',label='G160M')
    ax2.plot( data['date'][G160M_B], data['x_shift'][G160M_B],'g.')
    ax2.xaxis.set_ticklabels( ['' for item in ax2.xaxis.get_ticklabels()] )

    ax3 = fig.add_subplot(3,1,3)
    ax3.plot( data['date'][G140L_A], data['x_shift'][G140L_A],'y.',label='G140L')
    ax3.plot( data['date'][G140L_B], data['x_shift'][G140L_B],'y.')
    ax3.set_xlabel('DATE [MJD]', fontsize=20, fontweight='bold')

    ax.legend(shadow=True, numpoints=1, loc='upper left')
    fig.suptitle('FUV SHIFT_DISP - FPPOS OFFSET [A/B]', fontsize=20, fontweight='bold')
    ax.set_ylabel('SHIFT1[A/B] (pixels)', fontsize=10, fontweight='bold')

    for axis,index in zip([ax,ax2,ax3],[G130M,G160M,G140L]):
        #axis.set_ylim(-300,300)
        axis.set_xlim(data['date'].min(),data['date'].max()+50 )
        axis.set_ylabel('SHIFT_DISP - FPPOS OFFSET [A/B]', fontsize=7, fontweight='bold')
        axis.axhline(y=0,color='r')
        axis.axhline(y=285,color='k',lw=3,ls='--',zorder=1,label='Search Range')
        axis.axhline(y=-285,color='k',lw=3,ls='--',zorder=1)
        fit,ydata,parameters,err = fit_data(data['date'][index],data['x_shift'][index])
        axis.plot( ydata,fit,'k-',lw=3,label='%3.5fx'%(parameters[0]) )
        axis.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1, numpoints=1,shadow=True,prop={'size':10})

    remove_if_there(os.path.join(out_dir,'FUV_shifts.png'))
    fig.savefig(os.path.join(out_dir,'FUV_shifts.png'))
    plt.close(fig)
    os.chmod(os.path.join(out_dir,'FUV_shifts.png'),0o766)

    ##########

    fig = plt.figure(figsize=(16, 18))
    ax = fig.add_subplot(7, 1, 1)
    ax.plot(data['date'][G185M_A].data, data['x_shift'][G185M_A].data, 'bo', label='G185M')
    ax.plot(data['date'][G185M_B].data, data['x_shift'][G185M_B].data, 'bo', markeredgecolor='k')
    ax.plot(data['date'][G185M_C].data, data['x_shift'][G185M_C].data, 'bo', markeredgecolor='k')
    ax.axhline(y=0, color='red')

    #--second timeframe
    transition_fraction = (56500.0 - data['date'].min()) / \
        (data['date'].max() - data['date'].min())

    ax.axhline(y=58, xmin=0, xmax=transition_fraction, color='k',
                lw=3, ls='--', zorder=1, label='Search Range')
    ax.axhline(y=-58, xmin=0, xmax=transition_fraction,
                color='k', lw=3, ls='--', zorder=1)

    ax.axhline(y=90, xmin=transition_fraction, xmax=1,
                color='k', lw=3, ls='--', zorder=1)
    ax.axhline(y=-90, xmin=transition_fraction,
                xmax=1, color='k', lw=3, ls='--', zorder=1)
    #--

    sigma = data['x_shift'][G185M_A].std()

    ax.xaxis.set_ticklabels(['' for item in ax.xaxis.get_ticklabels()])

    ax2 = fig.add_subplot(7, 1, 2)
    ax2.plot(data['date'][G225M_A], data['x_shift'][G225M_A], 'ro', label='G225M')
    ax2.plot(data['date'][G225M_B], data['x_shift'][G225M_B], 'ro', markeredgecolor='k')
    ax2.plot(data['date'][G225M_C], data['x_shift'][G225M_C], 'ro', markeredgecolor='k')
    ax2.axhline(y=0, color='red')

    #--second timeframe
    transition_fraction = (56500.0 - data['date'].min()) / \
        (data['date'].max() - data['date'].min())

    ax2.axhline(y=58, xmin=0, xmax=transition_fraction, color='k', lw=3, ls='--', zorder=1, label='Search Range')
    ax2.axhline(y=-58, xmin=0, xmax=transition_fraction, color='k', lw=3, ls='--', zorder=1)

    ax2.axhline(y=90, xmin=transition_fraction, xmax=1,
                color='k', lw=3, ls='--', zorder=1)
    ax2.axhline(y=-90, xmin=transition_fraction,
                xmax=1, color='k', lw=3, ls='--', zorder=1)
    #--

    sigma = data['x_shift'][G225M_A].std()

    ax2.xaxis.set_ticklabels(['' for item in ax2.xaxis.get_ticklabels()])

    ax3 = fig.add_subplot(7, 1, 3)
    ax3.plot(data['date'][G285M_A], data['x_shift'][G285M_A], 'yo', label='G285M')
    ax3.plot(data['date'][G285M_B], data['x_shift']
             [G285M_B], 'yo', markeredgecolor='k')
    ax3.plot(data['date'][G285M_C], data['x_shift']
             [G285M_C], 'yo', markeredgecolor='k')
    ax3.axhline(y=0, color='red')
    ax3.axhline(y=58, color='k', lw=3, ls='--', zorder=1, label='Search Range')
    ax3.axhline(y=-58, color='k', lw=3, ls='--', zorder=1)

    sigma = data['x_shift'][G285M_A].std()

    ax3.xaxis.set_ticklabels(['' for item in ax3.xaxis.get_ticklabels()])

    ax4 = fig.add_subplot(7, 1, 4)
    ax4.plot(data['date'][G230L_A], data['x_shift'][G230L_A], 'go', label='G230L')
    ax4.plot(data['date'][G230L_B], data['x_shift']
             [G230L_B], 'go', markeredgecolor='k')
    ax4.plot(data['date'][G230L_C], data['x_shift']
             [G230L_C], 'go', markeredgecolor='k')

    ax4.axhline(y=0, color='red')

    #--second timeframe
    transition_fraction = (55535.0 - data['date'].min()) / \
        (data['date'].max() - data['date'].min())

    ax4.axhline(y=58, xmin=0, xmax=transition_fraction, color='k',
                lw=3, ls='--', zorder=1, label='Search Range')
    ax4.axhline(y=-58, xmin=0, xmax=transition_fraction,
                color='k', lw=3, ls='--', zorder=1)

    ax4.axhline(y=90, xmin=transition_fraction, xmax=1,
                color='k', lw=3, ls='--', zorder=1)
    ax4.axhline(y=-90, xmin=transition_fraction,
                xmax=1, color='k', lw=3, ls='--', zorder=1)
    #--
    ax4.xaxis.set_ticklabels(['' for item in ax3.xaxis.get_ticklabels()])
    sigma = data['x_shift'][G230L_A].std()

    ax.set_title('NUV SHIFT_DISP - FPPOS OFFSET [A/B/C]', fontsize=20, fontweight='bold')
    for axis, index in zip([ax, ax2, ax3, ax4], [G185M, G225M, G285M, G230L]):
        #axis.set_ylim(-110, 110)
        axis.set_xlim(data['date'].min(), data['date'].max() + 50)
        axis.set_ylabel('SHIFT_DISP - FPPOS OFFSET [A/B]', fontsize=7, fontweight='bold')
        fit, ydata, parameters, err = fit_data(
            data['date'][index], data['x_shift'][index])
        axis.plot(ydata, fit, 'k-', lw=3, label='%3.5fx' % (parameters[0]))
        axis.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1, numpoints=1, shadow=True, fontsize=12)

    ax4.set_xlabel('date',fontsize=20, fontweight='bold')

    ax = fig.add_subplot(7, 1, 5)
    ax.plot(data['date'][NUV], data['x_shift'][NUV], '.')
    fit, ydata, parameters, err = fit_data(
        data['date'][NUV], data['x_shift'][NUV])
    ax.plot(ydata, fit, 'k-', lw=3, label='%3.5fx' % (parameters[0]))
    ax.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1,numpoints=1, shadow=True)
    ax.set_ylabel('All NUV', fontsize=15, fontweight='bold')
    ax.xaxis.set_ticklabels(['' for item in ax.xaxis.get_ticklabels()])
    ax.set_xlim(data['date'].min(), data['date'].max() + 50)
    #ax.set_ylim(-110, 110)
    
    mirrora = np.where((data_acqs['opt_elem'] == 'MIRRORA')
                       & (data_acqs['x_shift'] > 0))[0]
    ax = fig.add_subplot(7, 1, 6)
    ax.plot(data_acqs['date'][mirrora], data_acqs['x_shift'][mirrora], '.')
    fit, ydata, parameters, err = fit_data(
        data_acqs['date'][mirrora], data_acqs['x_shift'][mirrora])
    ax.plot(ydata, fit, 'k-', lw=3, label='%3.5fx' % (parameters[0]))
    ax.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1,numpoints=1, shadow=True)
    ax.set_xlim(data_acqs['date'].min(), data_acqs['date'].max() + 50)
    ax.set_ylabel('MIRRORA', fontsize=15, fontweight='bold')
    #ax.set_ylim(460, 630)

    mirrorb = np.where((data_acqs['opt_elem'] == 'MIRRORB')
                       & (data_acqs['x_shift'] > 0))[0]
    ax = fig.add_subplot(7, 1, 7)
    ax.plot(data_acqs['date'][mirrorb], data_acqs['x_shift'][mirrorb], '.')
    fit, ydata, parameters, err = fit_data(
        data_acqs['date'][mirrorb], data_acqs['x_shift'][mirrorb])
    ax.plot(ydata, fit, 'k-', lw=3, label='%3.5fx' % (parameters[0]))
    ax.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1,numpoints=1, shadow=True)
    ax.set_xlim(data_acqs['date'].min(), data_acqs['date'].max() + 50)
    ax.set_ylabel('MIRRORB', fontsize=15, fontweight='bold')
    ax.set_xlabel('Date [MJD]', fontsize=20, fontweight='bold')
    #ax.set_ylim(260, 400)

    remove_if_there(os.path.join(out_dir, 'NUV_shifts.png'))
    fig.savefig(os.path.join(out_dir, 'NUV_shifts.png'),
                bbox_inches='tight',
                pad_inches=.5)
    plt.close(fig)
    os.chmod(os.path.join(out_dir, 'NUV_shifts.png'),0o766)

    ##############

    for elem in ['MIRRORA', 'MIRRORB']:
        mirror = np.where((data_acqs['opt_elem'] == elem)
                          & (data_acqs['x_shift'] > 0))[0]
        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(data_acqs['date'][mirror], data_acqs['x_shift'][mirror], '.')
        fit, ydata, parameters, err = fit_data(data_acqs['date'][mirror],
                                               data_acqs['x_shift'][mirror])
        ax.plot(ydata, fit, 'r-', lw=3, label='%3.5f +/- %3.5f' %
                (parameters[0], err))
        ax.legend(numpoints=1, shadow=True, loc='upper left')
        ax.set_xlim(data_acqs['date'].min(), data_acqs['date'].max() + 50)
        #ax.set_ylim(460, 630)
        remove_if_there(os.path.join(out_dir, '{}_shifts.png'.format(elem.upper())))
        fig.savefig(os.path.join(out_dir, '{}_shifts.png'.format(elem.upper())))
        plt.close(fig)
        os.chmod((os.path.join(out_dir, '{}_shifts.png'.format(elem.upper()))),0o766)


    for grating in list(set(data['opt_elem'])):
        fig = plt.figure()
        ax = fig.add_axes([.1, .1, .75, .8])
        ax.set_title(grating, fontsize=20, fontweight='bold')
        for cenwave in list(set(data['cenwave'])):
            index = np.where((data['opt_elem'] == grating) &
                             (data['cenwave'] == cenwave))[0]
            if not len(index):
                continue

            xdata = np.array(map(int, data['date'][index]))
            ydata = data['x_shift'][index]
            new_ydata = []
            new_xdata = []
            for day in range(xdata.min(), xdata.max() + 1):
                index = np.where(xdata == day)[0]
                #n_times = len(index)
                median = np.median(ydata[index])
                new_ydata.append(median)
                new_xdata.append(day)

            if cenwave < 1700:
                ms = 6
                ylim = (-140, 80)
            else:
                ms = 10
                ylim = (-80, 80)

            ax.plot(new_xdata, new_ydata, '.', ms=ms, alpha=.7, label='%d' %
                    (cenwave))

            plt.legend(numpoints=1, shadow=True, bbox_to_anchor=(1.05, 1),
                       loc='upper left', borderaxespad=0., prop={'size': 8})
            ax.set_xlim(data['date'].min(), data['date'].max() + 50)
            #ax.set_ylim(ylim[0], ylim[1])
        remove_if_there(os.path.join(out_dir, '%s_shifts_color.pdf' %
                    (grating)))
        fig.savefig(os.path.join(out_dir, '%s_shifts_color.pdf' %
                    (grating)))
        plt.close(fig)
        os.chmod(os.path.join(out_dir, '%s_shifts_color.pdf' %
                    (grating)), 0o766)

#----------------------------------------------------------

def make_plots_2(data, data_acqs, out_dir):
    """ Making the plots for the shift2 value
    """

    sorted_index = np.argsort(data['date'])
    data = data[sorted_index]

    for cenwave in set(data['cenwave']):
        cw_index = np.where(data['cenwave'] == cenwave)
        all_segments = set(data[cw_index]['segment'])
        n_seg = len(all_segments)

        fig = plt.figure()
        fig.suptitle('Shift2 vs Shift1 {}'.format(cenwave), fontsize=20, fontweight='bold')

        for i, segment in enumerate(all_segments):
            index = np.where( (data['segment'] == segment) &
                              (data['cenwave'] == cenwave) )

            ax = fig.add_subplot(n_seg, 1, i+1)
            ax.plot(data[index]['x_shift'], data[index]['y_shift'], 'o')
            ax.set_xlabel('x_shift', fontsize=20, fontweight='bold')
            ax.set_ylabel('y_shift', fontsize=15, fontweight='bold')
            #ax.set_ylabel('SHIFT2 vs SHIFT1 {}'.format(segment))
            #ax.set_ylim(-20, 20)
        remove_if_there(os.path.join(out_dir, 'shift_relation_{}.png'.format(cenwave)))
        fig.savefig(os.path.join(out_dir, 'shift_relation_{}.png'.format(cenwave)))
        plt.close(fig)
        os.chmod(os.path.join(out_dir, 'shift_relation_{}.png'.format(cenwave)), 0o766)


#----------------------------------------------------------

def fp_diff(data):
    index = np.where((data['detector'] == 'FUV'))[0]
    data = data[index]

    datasets = list(set(data['dataset']))
    datasets.sort()

    all_cenwaves = set(data['cenwave'])
    diff_dict = {}
    for cenwave in all_cenwaves:
        diff_dict[cenwave] = []

    ofile = open(os.path.join(out_dir, 'shift_data.txt'), 'w')
    for name in datasets:
        a_shift = None
        b_shift = None
        try:
            a_shift = data['x_shift'][np.where((data['dataset'] == name) &
                                               (data['segment'] == 'FUVA'))[0]][0]
            b_shift = data['x_shift'][np.where((data['dataset'] == name) &
                                               (data['segment'] == 'FUVB'))[0]][0]
        except IndexError:
            continue

        cenwave = data['cenwave'][np.where((data['dataset'] == name) &
                                           (data['segment'] == 'FUVA'))[0]][0]
        opt_elem = data['opt_elem'][np.where((data['dataset'] == name) &
                                             (data['segment'] == 'FUVA'))[0]][0]
        fppos = data['fppos'][np.where((data['dataset'] == name) &
                                       (data['segment'] == 'FUVA'))[0]][0]
        mjd = data['date'][np.where((data['dataset'] == name) &
                                   (data['segment'] == 'FUVA'))[0]][0]
        diff = a_shift - b_shift

        diff_dict[cenwave].append((mjd, diff))
        ofile.write('%5.5f  %s  %d  %d   %3.2f  %3.2f  \n' %
                    (mjd, opt_elem, cenwave, fppos, a_shift, b_shift))

    for cenwave in diff_dict:
        all_diff = [line[1] for line in diff_dict[cenwave]]
        all_mjd = [line[0] for line in diff_dict[cenwave]]

        if not len(all_diff):
            continue

        plt.figure(figsize=(8, 5))
        plt.plot(all_mjd, all_diff, 'o', label='%s' % (cenwave))
        plt.xlabel('MJD')
        plt.ylabel('SHIFT1 difference (pixels)')
        plt.title(cenwave, fontsize=20, fontweight='bold')
        plt.legend(shadow=True, numpoints=1, loc='upper left')
        remove_if_there(os.path.join(out_dir, 'difference_%s.pdf' % (cenwave)))
        plt.savefig(os.path.join(out_dir, 'difference_%s.pdf' % (cenwave)))
        plt.close()
        os.chmod(os.path.join(out_dir, 'difference_%s.pdf' % (cenwave)), 0o766)

#----------------------------------------------------------

def monitor():
    """Run the entire suite of monitoring
    """

    logger.info("STARTING MONITOR")

    settings = get_settings()

    webpage_dir = os.path.join(settings['webpage_location'], 'shifts')
    monitor_dir = os.path.join(settings['monitor_location'], 'Shifts')

    for place in [webpage_dir, monitor_dir]:
        if not os.path.exists(place):
            logger.debug("CREATING MONITOR LOCATION: {}".format(place))
            os.makedirs(place)

    #-- Create tables
    flash_data = make_shift_table(Lampflash)
    rawacq_data = make_shift_table(Rawacqs)
    
    #-- Make static plots.
    make_plots(flash_data, rawacq_data, monitor_dir)
    
    make_interactive_plots(flash_data, rawacq_data, monitor_dir, 'FUV')
    make_interactive_plots(flash_data, rawacq_data, monitor_dir, 'NUV')

    # make_plots_2(flash_data, rawacq_data, monitor_dir)
    # fp_diff(flash_data)

    for item in glob.glob(os.path.join(monitor_dir, '*.p??')):
        remove_if_there(os.path.join(webpage_dir, os.path.basename(item)))
        shutil.copy(item, webpage_dir)

    logger.info("FINISH MONITOR")

#----------------------------------------------------------