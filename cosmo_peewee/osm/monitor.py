""" Script to compile the spectrum shift data for COS FUV and NUV data.
"""

from __future__ import absolute_import

import glob
import logging
import os
import shutil
import sys

from astropy import table
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.time import Time
from bokeh.io import output_file, show, save
from bokeh.layouts import column, row
from bokeh.models import Range1d, HoverTool, BoxSelectTool, ColumnDataSource, \
                         OpenURL, TapTool, Div, Button, CustomJS
from bokeh.plotting import figure
from copy import deepcopy
from datetime import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
from scipy.stats import linregress
from time import gmtime, strftime, localtime

from ..database.models import get_database, get_settings
from ..database.models import Lampflash, Rawacqs, Files
from ..utils import remove_if_there

logger = logging.getLogger(__name__)
np.seterr(divide='ignore')
np.warnings.filterwarnings('ignore') 
# Sometime the np.mean gets empty slice, its okay, just want to silence output.

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

    lamptab = fits.getdata(os.path.join(os.environ['lref'], 
                           lamptab_name))

    if 'FPOFFSET' not in lamptab.names:
        return 0

    index = np.where((lamptab['segment'] == segment)
                     & (lamptab['opt_elem'] == opt_elem)
                     & (lamptab['cenwave'] == cenwave)
                     & (lamptab['fpoffset'] == fpoffset))[0]

    offset = lamptab['FP_PIXEL_SHIFT'][index][0]

    return offset

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

    # Open file
    file_path = os.path.join(filename.path, filename.filename)
    with fits.open(file_path) as hdu:
        # Set some dictionary values.
        out_info = {'filename': filename.filename,
                    'date': round(hdu[1].header['EXPSTART'],5),
                    'rootname': hdu[0].header['ROOTNAME'],
                    'proposid': hdu[0].header['PROPOSID'],
                    'detector': hdu[0].header['DETECTOR'],
                    'opt_elem': hdu[0].header['OPT_ELEM'],
                    'cenwave': hdu[0].header['CENWAVE'],
                    'fppos': hdu[0].header.get('FPPOS', None),
                    'filetype': hdu[0].header.get('FILETYPE', None)}

        # Get time, and then convert the format
        t = Time(out_info['date'], format='mjd')
        out_info['cal_date'] = t.iso

        # Open lampflash
        if '_lampflash.fits' in filename.filename:
            out_info['life_adj'] = hdu[0].header.get('LIFE_ADJ')
            out_info['segment'] = hdu[0].header['SEGMENT']
            # Get lamptab file
            out_info['lamptab'] = hdu[0].header['LAMPTAB'].split('$')[-1]

            # FPPOS 3 is the home frame, so put all FP's in home frame.
            fpoffset = out_info['fppos'] - 3

            for i, line in enumerate(hdu[1].data):
                # 'flash' counts the number of flashes in a lampflash
                # x_shift (dispersion axis) if the calculated shift for the monitor.
                out_info['flash'] = (i // 2) + 1
                out_info['x_shift'] = line['SHIFT_DISP'] \
                                      - fppos_shift(out_info['lamptab'],
                                                    line['segment'],
                                                    out_info['opt_elem'],
                                                    out_info['cenwave'],
                                                    fpoffset)

                out_info['y_shift'] = line['SHIFT_XDISP']
                out_info['found'] = line['SPEC_FOUND']
                out_info['segment'] = line['SEGMENT']

                # don't need too much precision here
                out_info['x_shift'] = round(out_info['x_shift'], 5)
                out_info['y_shift'] = round(out_info['y_shift'], 5)

                yield deepcopy(out_info)
        
        # Open rawacqs
        elif '_rawacq.fits' in filename.filename:
            #-- Technically it wasn't found.
            out_info['found'] = False
            out_info['fppos'] = -1
            out_info['flash'] = 1

            # Grab associated spt
            spt = fits.open(os.path.join(filename.path,
                                         filename.filename.replace('rawacq', 
                                                                   'spt')))
            
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


def fit_data(xdata, ydata):
    """ Fit a regression line to shift data points

    Parameters
    ----------
    xdata : astropy.table.column.Column
        A list of x values
    ydata : astropy.table.column.Column
        A list of y values

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

    # this is a crude implementation, but it lets me use the rest of the
    # plotting code as-is

    # I dont think this will work with python 3 :(
    # .dicts() returns the result objects as dictionaries. 
    for i, row in enumerate(db_table.select().dicts()):
        data.append(row.values())
        if not i:
            # get keys here because if you use ._meta.fields.keys() 
            # they will be out of order.
            keys = row.keys()    
    
    database.close()

    data = Table(rows=data, names=keys)

    return data

def make_panel(data, **kwargs):
    """Make a bokeh panel for figure.

    Parameters
    ----------
    data: Astropy.Table
        Astropy table of all metadata.
    **kwargs
        Keyword arguements for observing configuration and plotting

        Depending on the configuration and location of the panel in the plot
        you will want to pass specific parameters to the bokeh panel object.
        For instance, the top panel will contain the title and the different
        panels will have different colors etc. Easier to cut down on the 
        number of arguements.
    """

    # Define tools that each panel will possess.
    TOOLS ='box_zoom,box_select,pan,reset,tap'
    
    grating = kwargs.get('grating')
    height = kwargs.get('height')
    width = kwargs.get('width')
    detector = kwargs.get('detector')
    plt_color = kwargs.get('plt_color')
    top = kwargs.get('top', False)
    x_range = kwargs.get('x_range', False)
    acqs = kwargs.get('acqs', False)
    all_nuv = kwargs.get('all_nuv', False)

    # Build ColumnDataSource object
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

    # Parse detector for labeling
    if detector == 'FUV':
        # If top panel, add a title
        if top:
            panel = figure(width=width, height=height, 
                           x_range=(min(data['date']) - 10, 
                                    max(data['date']) + 10), 
                           title='FUV SHIFT1[A/B] as of {} EST'   
                                 .format(strftime("%m-%d-%Y %H:%M:%S", 
                                                  localtime())), 
                           tools=[TOOLS,hover])
            panel.title.text_font_size = '15pt'
        else:
            panel = figure(width=width, height=height, x_range=x_range, 
                           title=None, tools=[TOOLS, hover])
       
        # Label y and also draw max bounds from reference file. 
        panel.line(data['date'], np.zeros_like(data['date']) + 300, 
                   color='black', line_width=2, line_dash='dashed')
        panel.line(data['date'], np.zeros_like(data['date']) - 300, 
                   color='black', line_width=2, line_dash='dashed')
        panel.yaxis.axis_label = "Shift1[A/B] (Pixels)"
    
    elif detector == 'NUV':
        # If top panel, add a title
        if top:
            panel = figure(width=width, height=height, 
                           x_range=(min(data['date']) - 10, 
                                    max(data['date']) + 10), 
                           title='NUV SHIFT1[A/B/C] as of {} EST'
                                 .format(strftime("%m-%d-%Y %H:%M:%S", 
                                         localtime())), 
                           tools=[TOOLS, hover])
            panel.title.text_font_size = '15pt'            
        else:
            panel = figure(width=width, height=height, x_range=x_range, 
                           title=None, tools=[TOOLS, hover])
        panel.yaxis.axis_label = "Shift1[A/B/C] (Pixels)"

    # Make scatter plot of data
    if all_nuv:
        panel.circle('date', 'shift', legend='All NUV', size=4, source=source, 
                     color=plt_color, alpha=0.5)
    else:
        panel.circle('date', 'shift', legend=data['opt_elem'][grating][0], 
                     size=4, source=source, color=plt_color, alpha=0.5)

    # Provide URL and taptool and callback info.
    url = "http://archive.stsci.edu/proposal_search.php?id=@proposid&mission=hst"
    taptool = panel.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    
    return panel

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
    
    Returns
    -------
    p: bokeh plot
        Bokeh plot object.
    """
    
    logger.info("MAKING INTERACTIVE PLOT FOR {}".format(detector))
    settings = get_settings()
    # Sort by time
    sorted_index = np.argsort(data['date'])
    
    # Rearange Data
    data = data[sorted_index]

    if detector == 'FUV':
        # Sort all data by opt_elem
        G140L = np.where((data['opt_elem'] == 'G140L'))[0]
        G130M = np.where((data['opt_elem'] == 'G130M'))[0]
        G160M = np.where((data['opt_elem'] == 'G160M'))[0]

        # Find Unique Entries
        unique_data = table.unique(data, ['date','x_shift',
                                          'rootname','proposid'])
        unique_G140L = np.where((unique_data['opt_elem'] == 'G140L'))[0]
        unique_G130M = np.where((unique_data['opt_elem'] == 'G130M'))[0]
        unique_G160M = np.where((unique_data['opt_elem'] == 'G160M'))[0]

        # Begin Bokeh   
        TOOLS ='box_zoom,box_select,crosshair,pan,reset,tap'
        
        # Plot FUV Shifts
        outname = os.path.join(out_dir, 'FUV_shift_vs_time.html')
        remove_if_there(outname)
        output_file(outname)
        
        # Set panel size
        plt_hgt = 250
        plt_wth = 800

        # Create bokeh figure objects
        
        # Panel 1    
        s1 = make_panel(unique_data, grating=unique_G130M, height=plt_hgt,
                        width=plt_wth, detector='FUV', plt_color='blue', 
                        top=True)
        # Fit shift as a function of date and plot it...
        fit,ydata,parameters,err = fit_data(data['date'][G130M],
                                            data['x_shift'][G130M])
        s1.line(ydata, fit, color='black', line_width=2, 
                legend=str(parameters[0]))
        
        # Panel 2
        s2 = make_panel(unique_data, grating=unique_G160M, height=plt_hgt,
                        width=plt_wth, detector='FUV', plt_color='green', 
                        x_range=s1.x_range)
        # Fit shift as a function of date and plot it...
        fit,ydata,parameters,err = fit_data(data['date'][G160M],
                                            data['x_shift'][G160M])
        s2.line(ydata, fit, color='black', line_width=2, 
                legend=str(parameters[0]))
        
        # Panel 3
        s3 = make_panel(unique_data, grating=unique_G140L, height=plt_hgt, 
                        width=plt_wth, detector='FUV', plt_color='yellow', 
                        x_range=s1.x_range)
        
        # Fit shift as a function of date and plot it...
        fit,ydata,parameters,err = fit_data(data['date'][G140L],
                                            data['x_shift'][G140L])
        s3.line(ydata, fit, color='black', line_width=2, 
                legend=str(parameters[0]))
        
        s3.xaxis.axis_label = "Time (MJD)"
        
        p = column(s1, s2, s3)

        save(p, filename=outname)

        return p
        
    # NUV Plots
    if detector == 'NUV':
        
        # Sort by grating.
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


        # Sort unique entries by grating.
        unique_data = table.unique(data, 
                                   ['date','x_shift','rootname','proposid'])
        unique_acqs = table.unique(data_acqs, 
                                   ['date','x_shift','rootname','proposid'])

        unique_G230L = np.where((unique_data['opt_elem'] == 'G230L'))[0]
        unique_G225M = np.where((unique_data['opt_elem'] == 'G225M'))[0]
        unique_G285M = np.where((unique_data['opt_elem'] == 'G285M'))[0]
        unique_G185M = np.where((unique_data['opt_elem'] == 'G185M'))[0]

        unique_G230L_A = np.where((unique_data['opt_elem'] == 'G230L') &
                                  (unique_data['segment'] == 'NUVA'))[0]
        unique_G230L_B = np.where((unique_data['opt_elem'] == 'G230L') &
                                  (unique_data['segment'] == 'NUVB'))[0]
        unique_G230L_C = np.where((unique_data['opt_elem'] == 'G230L') &
                                  (unique_data['segment'] == 'NUVC'))[0]


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
        
        wcptab_path = os.path.join(settings['lref'], '03p1706jl_wcp.fits')
        wcptab = fits.getdata(wcptab_path)
        wcptab_data = [row for row in wcptab]
        names = ['OPT_ELEM', 'XC_RANGE', 'RESWIDTH', 'MAX_TIME_DIFF',
                 'STEPSIZE', 'XD_RANGE', 'BOX', 'SEARCH_OFFSET']
        wcptab_table = Table(rows=wcptab_data, names=names)

        # Bokeh
        TOOLS ='box_zoom,pan,reset,hover,tap'

        # Bokeh panel sizes.
        plt_hgt = 250
        plt_wth = 800

        # Set outname and create file.
        outname = os.path.join(out_dir, 'NUV_shift_vs_time.html')
        remove_if_there(outname)
        output_file(outname)
        
        # G230L search range was updated earlier than the other observing modes.
        transition_date = 56500.0
        transition_date_G230L = 55535.0        

        # Because of the complexity of the different transition dates, 
        # plotting with bokeh and acq figures... I've made code blocks 
        # for each panel. 
        
        # Panel 1
        # Create bokeh figure.
        s1 = make_panel(unique_data, grating=unique_G185M, height=plt_hgt, 
                        width=plt_wth, detector='NUV', plt_color='blue', 
                        top=True)
        # Fit Data
        fit,ydata,parameters,err = fit_data(data['date'][G185M],
                                            data['x_shift'][G185M])
        s1.line(ydata, fit, color='black', 
                line_width=2, legend=str(parameters[0]))

        # Find transition regions.
        before_data = np.where(data['date'][G185M] <= transition_date)
        after_data = np.where(data['date'][G185M] >= transition_date)
        
        wcptab_row = np.where(wcptab_table['OPT_ELEM']=='G185M')
        xc_range = wcptab_table[wcptab_row]['XC_RANGE']
        srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']

        # First transitions
        s1.line(data['date'][G185M][before_data], 
                np.zeros_like(data['date'][G185M][before_data]) 
                + 58, color='black', line_width=2, line_dash='dashed')
        s1.line(data['date'][G185M][before_data], 
                np.zeros_like(data['date'][G185M][before_data]) 
                - 58, color='black', line_width=2, line_dash='dashed')

        # Second
        s1.line(data['date'][G185M][after_data], 
                np.zeros_like(data['date'][G185M][after_data]) 
                + (xc_range + srh_offset), 
                color='black', line_width=2, line_dash='dashed')
        s1.line(data['date'][G185M][after_data], 
                np.zeros_like(data['date'][G185M][after_data]) 
                + ((xc_range *-1) + srh_offset), 
                color='black', line_width=2, line_dash='dashed')
        
        s1.line(data['date'][G185M][before_data], 
                np.zeros_like(data['date'][G185M][before_data]), 
                color='red', line_width=2)
        s1.line(data['date'][G185M][after_data], 
                np.zeros_like(data['date'][G185M][after_data]) 
                + srh_offset, color='red', line_width=2)

        # Panel 2
        wcptab_row = np.where(wcptab_table['OPT_ELEM']=='G225M')
        xc_range = wcptab_table[wcptab_row]['XC_RANGE']
        srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']

        s2 = make_panel(unique_data, grating=unique_G225M, height=plt_hgt, 
                        width=plt_wth, detector='NUV', plt_color='red', 
                        x_range=s1.x_range)

        fit,ydata,parameters,err = fit_data(data['date'][G225M],
                                            data['x_shift'][G225M])
        s2.line(ydata, fit, color='black', 
                line_width=2, legend=str(parameters[0]))
        before_data = np.where(data['date'][G225M] <= transition_date)
        after_data = np.where(data['date'][G225M] >= transition_date)
        
        s2.line(data['date'][G225M][before_data], 
                np.zeros_like(data['date'][G225M][before_data]) 
                + 58, color='black', line_width=2, line_dash='dashed')
        s2.line(data['date'][G225M][before_data], 
                np.zeros_like(data['date'][G225M][before_data]) 
                - 58, color='black', line_width=2, line_dash='dashed')

        s2.line(data['date'][G225M][after_data], 
                np.zeros_like(data['date'][G225M][after_data]) 
                + (xc_range + srh_offset), 
                color='black', line_width=2, line_dash='dashed')
        s2.line(data['date'][G225M][after_data], 
                np.zeros_like(data['date'][G225M][after_data]) 
                + ((xc_range *-1) + srh_offset), 
                color='black', line_width=2, line_dash='dashed')
        
        s2.line(data['date'][G225M][before_data], 
                np.zeros_like(data['date'][G225M][before_data]), 
                color='red', line_width=2)
        s2.line(data['date'][G225M][after_data], 
                np.zeros_like(data['date'][G225M][after_data]) 
                + srh_offset, color='red', line_width=2)

        # Panel 3
        wcptab_row = np.where(wcptab_table['OPT_ELEM']=='G285M')
        xc_range = wcptab_table[wcptab_row]['XC_RANGE']
        srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']
        
        s3 = make_panel(unique_data, grating=unique_G285M, height=plt_hgt, 
                        width=plt_wth, detector='NUV', plt_color='yellow', 
                        x_range=s1.x_range)
        
        fit,ydata,parameters,err = fit_data(data['date'][G285M],
                                            data['x_shift'][G285M])
        s3.line(ydata, fit, color='black', 
                line_width=2, legend=str(parameters[0]))
        
        before_data = np.where(data['date'][G285M] <= transition_date)
        after_data = np.where(data['date'][G285M] >= transition_date)
        s3.line(data['date'][G285M][before_data], 
                np.zeros_like(data['date'][G285M][before_data]) 
                + 58, color='black', line_width=2, line_dash='dashed')
        s3.line(data['date'][G285M][before_data], 
                np.zeros_like(data['date'][G285M][before_data]) 
                - 58, color='black', line_width=2, line_dash='dashed')

        s3.line(data['date'][G285M][after_data], 
                np.zeros_like(data['date'][G285M][after_data]) 
                + (xc_range + srh_offset), 
                color='black', line_width=2, line_dash='dashed')
        s3.line(data['date'][G285M][after_data], 
                np.zeros_like(data['date'][G285M][after_data]) 
                + ((xc_range *-1) + srh_offset), 
                color='black', line_width=2, line_dash='dashed')
        
        s3.line(data['date'][G285M][before_data], 
                np.zeros_like(data['date'][G285M][before_data]), 
                color='red', line_width=2)
        s3.line(data['date'][G285M][after_data], 
                np.zeros_like(data['date'][G285M][after_data]) 
                + srh_offset, color='red', line_width=2)
        
        # Panel 4
        wcptab_row = np.where(wcptab_table['OPT_ELEM']=='G230L')
        xc_range = wcptab_table[wcptab_row]['XC_RANGE']
        srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']

        s4 = make_panel(unique_data, grating=unique_G230L, height=plt_hgt, 
                        width=plt_wth, detector='NUV', plt_color='green', 
                        x_range=s1.x_range)

        fit,ydata,parameters,err = fit_data(data['date'][G230L],
                                            data['x_shift'][G230L])
        s4.line(ydata, fit, color='black', 
                line_width=2, legend=str(parameters[0]))

        before_data = np.where(data['date'][G230L] <= transition_date_G230L)
        after_data = np.where(data['date'][G230L] >= transition_date_G230L)
        s4.line(data['date'][G230L][before_data], 
                np.zeros_like(data['date'][G230L][before_data]) 
                + 58, color='black', line_width=2, line_dash='dashed')
        s4.line(data['date'][G230L][before_data], 
                np.zeros_like(data['date'][G230L][before_data]) 
                - 58, color='black', line_width=2, line_dash='dashed')

        s4.line(data['date'][G230L][after_data], 
                np.zeros_like(data['date'][G230L][after_data]) 
                + (xc_range + srh_offset), 
                color='black', line_width=2, line_dash='dashed')
        s4.line(data['date'][G230L][after_data], 
                np.zeros_like(data['date'][G230L][after_data]) 
                + ((xc_range *-1) + srh_offset), 
                color='black', line_width=2, line_dash='dashed')
        
        s4.line(data['date'][G230L][before_data], 
                np.zeros_like(data['date'][G230L][before_data]), 
                color='red', line_width=2)
        s4.line(data['date'][G230L][after_data], 
                np.zeros_like(data['date'][G230L][after_data]) 
                + srh_offset, color='red', line_width=2)
        
        # Panel 5
        s5 = make_panel(unique_data, grating=unique_NUV, height=plt_hgt, 
                        width=plt_wth, detector='NUV', plt_color='firebrick', 
                        x_range=s1.x_range, all_nuv=True)
        fit,ydata,parameters,err = fit_data(data['date'][NUV],
                                            data['x_shift'][NUV])
        s5.line(ydata, fit, color='black', 
                line_width=2, legend=str(parameters[0]))
        
        # Panel 6
        s6 = make_panel(unique_acqs, grating=unique_mirrora, height=plt_hgt, 
                        width=plt_wth, detector='NUV', plt_color='firebrick', 
                        x_range=s1.x_range, acqs=True)
        fit,ydata,parameters,err = fit_data(data_acqs['date'][mirrora],
                                            data_acqs['x_shift'][mirrora])
        s6.line(ydata, fit, color='black', 
                line_width=2, legend=str(parameters[0]))
        
        # Panel 7
        s7 = make_panel(unique_acqs, grating=unique_mirrorb, height=plt_hgt, 
                        width=plt_wth, detector='NUV', plt_color='firebrick', 
                        x_range=s1.x_range, acqs=True)
        fit,ydata,parameters,err = fit_data(data_acqs['date'][mirrorb],
                                            data_acqs['x_shift'][mirrorb])
        s7.line(ydata, fit, color='black', 
                line_width=2, legend=str(parameters[0]))
        s7.xaxis.axis_label = "Date (MJD)"
        
        # Format into single column.
        p = column(s1, s2, s3, s4, s5, s6, s7)
        save(p, filename=outname)

        return p        

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
    
    Returns
    -------
    None
    """
    
    logger.info("MAKING STATIC PLOTS")

    settings = get_settings()

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

    # Set plotting params.

    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(3,1,1)

    # Begin plotting FUV
    
    ax.plot(data['date'][G130M_A], data['x_shift'][G130M_A],
            'b+',label='G130M FUVA')
    ax.plot(data['date'][G130M_B], data['x_shift'][G130M_B],
            'rx',label='G130M FUVB')
    ax.xaxis.set_ticklabels( ['' for item in ax.xaxis.get_ticklabels()])

    ax2 = fig.add_subplot(3,1,2)
    ax2.plot(data['date'][G160M_A], data['x_shift'][G160M_A],
             'b+',label='G160M FUVA')
    ax2.plot(data['date'][G160M_B], data['x_shift'][G160M_B],
             'rx',label='G160M FUVB')
    ax2.xaxis.set_ticklabels( ['' for item in ax2.xaxis.get_ticklabels()])

    ax3 = fig.add_subplot(3,1,3)
    ax3.plot(data['date'][G140L_A], data['x_shift'][G140L_A],
             'b+',label='G140L FUVA')
    ax3.plot(data['date'][G140L_B], data['x_shift'][G140L_B],
             'rx',label='G140L FUVB')
    ax3.set_xlabel('DATE [MJD]', fontsize=20, fontweight='bold')

    ax.legend(shadow=True, numpoints=1, loc='upper left')
    fig.suptitle('FUV[A/B] SHIFT1', fontsize=20, fontweight='bold')
    ax.set_ylabel('SHIFT1[A/B] (pixels)', fontsize=10, fontweight='bold')

    # Change wcptab to fuv_wcptab in future to avoid var name conflict
    # with NUV.
    wcptab = fits.getdata(os.path.join(settings['lref'], '14o20140l_wcp.fits'))
    wcptab_data = [row for row in wcptab]
    names = ['OPT_ELEM', 'XC_RANGE', 'RESWIDTH', 'MAX_TIME_DIFF', 'STEPSIZE', 
             'XD_RANGE', 'BOX', 'SEARCH_OFFSET']
    wcptab_table = Table(rows=wcptab_data, names=names)

    for axis,index,grat in zip([ax,ax2,ax3],
                               [G130M,G160M,G140L],
                               ['G130M','G160M','G140L']):
        #axis.set_ylim(-300,300)
        axis.set_xlim(data['date'].min(),data['date'].max()+50 )
        axis.set_ylabel('SHIFT1 [A/B]', fontsize=15, fontweight='bold')
        
        wcptab_row = np.where(wcptab_table['OPT_ELEM']==grat)
        xc_range = wcptab_table[wcptab_row]['XC_RANGE']
        # srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']
        srh_offset = 0
        axis.axhline(y=srh_offset,color='r')
        axis.axhline(y=0,color='r')
        axis.axhline(y=xc_range + abs(srh_offset),color='k',
                     lw=3,ls='--',zorder=1,label='Search Range')
        axis.axhline(y=(xc_range*-1) - abs(srh_offset),color='k',
                     lw=3,ls='--',zorder=1)
        
        fit,ydata,parameters,err = fit_data(data['date'][index],
                                            data['x_shift'][index])
        axis.plot(ydata,fit,'k-',lw=3,label='%3.5fx'%(parameters[0]))
        axis.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1, 
                    numpoints=1, shadow=True,prop={'size':10})

    remove_if_there(os.path.join(out_dir,'FUV_new_ds_shifts.png'))
    fig.savefig(os.path.join(out_dir,'FUV_new_ds_shifts.png'))
    plt.close(fig)
    os.chmod(os.path.join(out_dir,'FUV_new_ds_shifts.png'),0o766)

    # Make table for current reference file
    wcptab = fits.getdata(os.path.join(settings['lref'], 
                                       '03p1706jl_wcp.fits'))
    wcptab_data = [row for row in wcptab]
    names = ['OPT_ELEM', 'XC_RANGE', 'RESWIDTH', 'MAX_TIME_DIFF', 'STEPSIZE', 
             'XD_RANGE', 'BOX', 'SEARCH_OFFSET']
    wcptab_table = Table(rows=wcptab_data, names=names)

    # Begin Plotting NUV
    # G185M
    fig = plt.figure(figsize=(16, 18))
    ax = fig.add_subplot(7, 1, 1)
    ax.plot(data['date'][G185M_A], data['x_shift'][G185M_A], 
            'rx', label='G185M NUVA')
    ax.plot(data['date'][G185M_B], data['x_shift'][G185M_B], 
            'go', mfc='none', label='G185M NUVB')
    ax.plot(data['date'][G185M_C], data['x_shift'][G185M_C], 
            'b+', label='G185M NUVC')
    
    # The WCPTAB has different search offsets and search ranges over time.
    # We will hardcode the old transitions because keeping track of the 
    # reference files is difficult.

    # First time frame
    transition_fraction = (56500.0 - data['date'].min()) / \
        (data['date'].max() - data['date'].min())

    ax.axhline(y=58, xmin=0, xmax=transition_fraction, color='k',
                lw=3, ls='--', zorder=1, label='Search Range')
    ax.axhline(y=-58, xmin=0, xmax=transition_fraction,
                color='k', lw=3, ls='--', zorder=1)
    ax.axhline(y=0, xmin=0, xmax=transition_fraction,
                color='red')
    
    #-- Second time frame
    wcptab_row = np.where(wcptab_table['OPT_ELEM']=='G185M')
    xc_range = wcptab_table[wcptab_row]['XC_RANGE']
    srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']
    
    ax.axhline(y=xc_range + abs(srh_offset), xmin=transition_fraction, xmax=1,
                color='k', lw=3, ls='--', zorder=1)
    ax.axhline(y=(xc_range*-1) - abs(srh_offset), xmin=transition_fraction,
                xmax=1, color='k', lw=3, ls='--', zorder=1)
    ax.axhline(y=srh_offset, xmin=transition_fraction, xmax=1,
                color='red')

    sigma = data['x_shift'][G185M_A].std()

    ax.xaxis.set_ticklabels(['' for item in ax.xaxis.get_ticklabels()])
    
    # G225M
    ax2 = fig.add_subplot(7, 1, 2)
    ax2.plot(data['date'][G225M_A], data['x_shift'][G225M_A], 
             'rx', label='G225M NUVA')
    ax2.plot(data['date'][G225M_B], data['x_shift'][G225M_B], 
             'go', mfc='none', label='G225M NUVB')
    ax2.plot(data['date'][G225M_C], data['x_shift'][G225M_C], 
             'b+', label='G225M NUVC')

    # First time frame
    transition_fraction = (56500.0 - data['date'].min()) / \
        (data['date'].max() - data['date'].min())

    ax2.axhline(y=58, xmin=0, xmax=transition_fraction, color='k', 
                lw=3, ls='--', zorder=1, label='Search Range')
    ax2.axhline(y=-58, xmin=0, xmax=transition_fraction, color='k', 
                lw=3, ls='--', zorder=1)

    #-- Second time frame
    wcptab_row = np.where(wcptab_table['OPT_ELEM']=='G225M')
    xc_range = wcptab_table[wcptab_row]['XC_RANGE']
    srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']
    
    ax2.axhline(y=xc_range + abs(srh_offset), xmin=transition_fraction, xmax=1,
                color='k', lw=3, ls='--', zorder=1)
    ax2.axhline(y=(xc_range*-1) - abs(srh_offset), xmin=transition_fraction,
                xmax=1, color='k', lw=3, ls='--', zorder=1)
    ax2.axhline(y=0, xmin=0, xmax=transition_fraction, color='red')
    ax2.axhline(y=srh_offset, xmin=transition_fraction, xmax=1,
                color='red')
    sigma = data['x_shift'][G225M_A].std()

    ax2.xaxis.set_ticklabels(['' for item in ax2.xaxis.get_ticklabels()])
    
    # G285M 
    ax3 = fig.add_subplot(7, 1, 3)
    ax3.plot(data['date'][G285M_A], data['x_shift'][G285M_A], 
             'rx', label='G285M NUVA')
    ax3.plot(data['date'][G285M_B], data['x_shift'][G285M_B], 
             'go', mfc='none', label='G285M NUVB')
    ax3.plot(data['date'][G285M_C], data['x_shift'][G285M_C], 
             'b+', label='G285M NUVC')
    ax3.axhline(y=0, color='red')
    
    # First time frame
    transition_fraction = (55535.0 - data['date'].min()) / \
        (data['date'].max() - data['date'].min())

    ax3.axhline(y=58, xmin=0, xmax=transition_fraction, color='k', 
                lw=3, ls='--', zorder=1, label='Search Range')
    ax3.axhline(y=-58, xmin=0, xmax=transition_fraction, color='k', 
                lw=3, ls='--', zorder=1)

    #-- Second time frame
    wcptab_row = np.where(wcptab_table['OPT_ELEM']=='G285M')
    xc_range = wcptab_table[wcptab_row]['XC_RANGE']
    srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']
    
    ax3.axhline(y=xc_range + srh_offset, xmin=transition_fraction, xmax=1,
                color='k', lw=3, ls='--', zorder=1)
    ax3.axhline(y=(xc_range*-1) + srh_offset, xmin=transition_fraction,
                xmax=1, color='k', lw=3, ls='--', zorder=1)
    ax3.axhline(y=0, xmin=0, xmax=transition_fraction, color='red')
    ax3.axhline(y=srh_offset, xmin=transition_fraction, xmax=1,
                color='red')

    sigma = data['x_shift'][G285M_A].std()

    ax3.xaxis.set_ticklabels(['' for item in ax3.xaxis.get_ticklabels()])

    # G230L
    ax4 = fig.add_subplot(7, 1, 4)
    ax4.plot(data['date'][G230L_A], data['x_shift'][G230L_A], 
             'rx', label='G230L NUVA')
    ax4.plot(data['date'][G230L_B], data['x_shift'][G230L_B], 
             'go', mfc='none', label='G230L NUVB')
    ax4.plot(data['date'][G230L_C], data['x_shift'][G230L_C], 
             'b+', label='G230L NUVC')

    # First time frame
    transition_fraction = (55535.0 - data['date'].min()) / \
        (data['date'].max() - data['date'].min())

    ax4.axhline(y=58, xmin=0, xmax=transition_fraction, color='k', 
                lw=3, ls='--', zorder=1, label='Search Range')
    ax4.axhline(y=-58, xmin=0, xmax=transition_fraction, color='k', 
                lw=3, ls='--', zorder=1)

    #-- Second time frame
    wcptab_row = np.where(wcptab_table['OPT_ELEM']=='G230L')
    xc_range = wcptab_table[wcptab_row]['XC_RANGE']
    srh_offset = wcptab_table[wcptab_row]['SEARCH_OFFSET']
    
    ax4.axhline(y=xc_range + srh_offset, xmin=transition_fraction, xmax=1,
                color='k', lw=3, ls='--', zorder=1)
    ax4.axhline(y=(xc_range*-1) + srh_offset, xmin=transition_fraction,
                xmax=1, color='k', lw=3, ls='--', zorder=1)
    ax4.axhline(y=0, xmin=0, xmax=transition_fraction, color='red')
    ax4.axhline(y=srh_offset, xmin=transition_fraction, xmax=1,
                color='red')

    ax4.xaxis.set_ticklabels(['' for item in ax3.xaxis.get_ticklabels()])
    sigma = data['x_shift'][G230L_A].std()

    ax.set_title('NUV[A/B/C] SHIFT1', fontsize=20, fontweight='bold')
    for axis, index in zip([ax, ax2, ax3, ax4], [G185M, G225M, G285M, G230L]):
        #axis.set_ylim(-110, 110)
        axis.set_xlim(data['date'].min(), data['date'].max() + 50)
        axis.set_ylabel('SHIFT1 [A/B/C]', fontsize=15, fontweight='bold')
        fit, ydata, parameters, err = fit_data(
            data['date'][index], data['x_shift'][index])
        axis.plot(ydata, fit, 'k-', lw=3, label='%3.5fx' % (parameters[0]))
        axis.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1, numpoints=1, 
                    shadow=True, fontsize=12)

    ax4.set_xlabel('date',fontsize=20, fontweight='bold')

    ax = fig.add_subplot(7, 1, 5)
    ax.plot(data['date'][NUV], data['x_shift'][NUV], '.')
    fit, ydata, parameters, err = fit_data(
        data['date'][NUV], data['x_shift'][NUV])
    ax.plot(ydata, fit, 'k-', lw=3, label='%3.5fx' % (parameters[0]))
    ax.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1,
              numpoints=1, shadow=True)
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
    ax.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1,
              numpoints=1, shadow=True)
    ax.set_xlim(data_acqs['date'].min(), data_acqs['date'].max() + 50)
    ax.set_ylabel('MIRRORA', fontsize=15, fontweight='bold')
    ax.set_xticks([])
    #ax.set_ylim(460, 630)

    mirrorb = np.where((data_acqs['opt_elem'] == 'MIRRORB')
                       & (data_acqs['x_shift'] > 0))[0]
    ax = fig.add_subplot(7, 1, 7)
    ax.plot(data_acqs['date'][mirrorb], data_acqs['x_shift'][mirrorb], '.')
    fit, ydata, parameters, err = fit_data(
        data_acqs['date'][mirrorb], data_acqs['x_shift'][mirrorb])
    ax.plot(ydata, fit, 'k-', lw=3, label='%3.5fx' % (parameters[0]))
    ax.legend(bbox_to_anchor=(1,1), loc='upper left', ncol=1,
              numpoints=1, shadow=True)
    ax.set_xlim(data_acqs['date'].min(), data_acqs['date'].max() + 50)
    ax.set_ylabel('MIRRORB', fontsize=15, fontweight='bold')
    ax.set_xlabel('Date [MJD]', fontsize=20, fontweight='bold')
    #ax.set_ylim(260, 400)

    remove_if_there(os.path.join(out_dir, 'NUV_new_ds_shifts.png'))
    fig.savefig(os.path.join(out_dir, 'NUV_new_ds_shifts.png'),
                bbox_inches='tight',
                pad_inches=.5)
    plt.close(fig)
    os.chmod(os.path.join(out_dir, 'NUV_new_ds_shifts.png'),0o766)

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
        remove_if_there(os.path.join(out_dir, 
                                     '{}_shifts.png'.format(elem.upper())))
        fig.savefig(os.path.join(out_dir, 
                                 '{}_shifts.png'.format(elem.upper())))
        plt.close(fig)
        os.chmod((os.path.join(out_dir, 
                               '{}_shifts.png'.format(elem.upper()))),0o766)

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


def plot(axrow, detector, cenwave, fuva_data=None, fuvb_data=None, nuva_data=None, nuvb_data=None, nuvc_data=None):
    all_fppos = ['fp1','fp2','fp3','fp4']
    
    if detector == 'FUV':
    
        fp_1a = np.where(fuva_data['fppos']==1)
        fp_2a = np.where(fuva_data['fppos']==2)
        fp_3a = np.where(fuva_data['fppos']==3)
        fp_4a = np.where(fuva_data['fppos']==4)

        for fp, marker, fp_label in zip([fp_1a, fp_2a, fp_3a, fp_4a], 
                                        ['*','+','x','v'], 
                                        all_fppos): 
            axrow[0].plot(fuva_data['date'][fp], fuva_data['x_shift'][fp],
                          'r{}'.format(marker),
                          label='FUVA | {} | {}'.format(cenwave, fp_label))
        
        axrow[0].set_xlabel('Date [MJD]', fontsize=20, fontweight='bold')
        axrow[0].set_ylabel('SHIFT1A [Pixels]', fontsize=20, 
                            fontweight='bold')
        axrow[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.01),
                        ncol=4, fancybox=True, shadow=True, fontsize=15)

        fp_1b = np.where(fuvb_data['fppos']==1)
        fp_2b = np.where(fuvb_data['fppos']==2)
        fp_3b = np.where(fuvb_data['fppos']==3)
        fp_4b = np.where(fuvb_data['fppos']==4)

        for fp, marker, fp_label in zip([fp_1b, fp_2b, fp_3b, fp_4b], 
                                        ['*','+','x','v'], 
                                        all_fppos): 
            axrow[1].plot(fuva_data['date'][fp], fuva_data['x_shift'][fp],
                          'b{}'.format(marker),
                          label='FUVA | {} | {}'.format(cenwave, fp_label))
        
        axrow[1].set_xlabel('Date [MJD]', fontsize=20, fontweight='bold')
        axrow[1].set_ylabel('SHIFT1A [Pixels]', fontsize=20, 
                            fontweight='bold')
        axrow[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.01),
                        ncol=4, fancybox=True, shadow=True, fontsize=15)

    elif detector == 'NUV':
        fp_1a = np.where(nuva_data['fppos']==1)
        fp_2a = np.where(nuva_data['fppos']==2)
        fp_3a = np.where(nuva_data['fppos']==3)
        fp_4a = np.where(nuva_data['fppos']==4)

        for fp, marker, fp_label in zip([fp_1a, fp_2a, fp_3a, fp_4a], 
                                        ['*','+','x','v'], 
                                        all_fppos): 
            axrow[0].plot(nuva_data['date'][fp], nuva_data['x_shift'][fp],
                          'r{}'.format(marker),
                          label='NUVA | {} | {}'.format(cenwave, fp_label))
        
        axrow[0].set_xlabel('Date [MJD]', fontsize=20, fontweight='bold')
        axrow[0].set_ylabel('SHIFT1A [Pixels]', fontsize=20, fontweight='bold')
        axrow[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.01),
                        ncol=4, fancybox=True, shadow=True, fontsize=15)
        
        fp_1b = np.where(nuvb_data['fppos']==1)
        fp_2b = np.where(nuvb_data['fppos']==2)
        fp_3b = np.where(nuvb_data['fppos']==3)
        fp_4b = np.where(nuvb_data['fppos']==4)

        for fp, marker, fp_label in zip([fp_1b, fp_2b, fp_3b, fp_4b], 
                                        ['*','+','x','v'], 
                                        all_fppos): 
            axrow[1].plot(nuvb_data['date'][fp], nuvb_data['x_shift'][fp],
                          'r{}'.format(marker),
                          label='NUVB | {} | {}'.format(cenwave, fp_label))

        axrow[1].set_xlabel('Date [MJD]', fontsize=20, fontweight='bold')
        axrow[1].set_ylabel('SHIFT1A [Pixels]', fontsize=20, fontweight='bold')
        axrow[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.01),
                        ncol=4, fancybox=True, shadow=True, fontsize=15)
        
        fp_1c = np.where(nuvc_data['fppos']==1)
        fp_2c = np.where(nuvc_data['fppos']==2)
        fp_3c = np.where(nuvc_data['fppos']==3)
        fp_4c = np.where(nuvc_data['fppos']==4)
        
        for fp, marker, fp_label in zip([fp_1c, fp_2c, fp_3c, fp_4c], 
                                        ['*','+','x','v'], 
                                        all_fppos): 
            axrow[2].plot(nuvc_data['date'][fp], nuvc_data['x_shift'][fp],
                          'g{}'.format(marker),
                          label='NUVC | {} | {}'.format(cenwave, fp_label))

        axrow[2].set_xlabel('Date [MJD]', fontsize=20, fontweight='bold')
        axrow[2].set_ylabel('SHIFT1A [Pixels]', fontsize=20, fontweight='bold')
        axrow[2].legend(loc='upper center', bbox_to_anchor=(0.5, 1.01),
                        ncol=4, fancybox=True, shadow=True, fontsize=15)


def make_plots_per_fppos(data, out_dir):
    """Make plot per time as with different FP-POS
    """
    sorted_index = np.argsort(data['date'])
    data = data[sorted_index]

    plt.rc('font', weight='bold')
    plt.rc('xtick.major', size=5, pad=7)
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15) 

    gratings = set(data['opt_elem'])
        
    for grating in gratings:
        index = np.where(data['opt_elem'] == grating)
        cenwaves = set(data[index]['cenwave'])
        segments = set(data[index]['segment'])    
        detector = list(set(data[index]['detector']))[0]

        n_cen = len(cenwaves)
        n_seg = len(segments)
        
        fig, axes = plt.subplots(n_cen, n_seg, figsize=(40,30))

        for ax, cenwave in zip(axes, cenwaves):
            if detector == 'FUV':
                fuva_data = np.where((data['segment'] == 'FUVA')
                                   & (data['cenwave'] == cenwave))

                fuvb_data = np.where((data['segment'] == 'FUVB')
                                   & (data['cenwave'] == cenwave))
                plot(ax, detector, cenwave, fuva_data=data[fuva_data], 
                     fuvb_data=data[fuvb_data])
            elif detector == 'NUV':
                nuva_data = np.where((data['segment'] == 'NUVA') 
                                   & (data['cenwave'] == cenwave))
                nuvb_data = np.where((data['segment'] == 'NUVB')
                                   & (data['cenwave'] == cenwave))
                nuvc_data = np.where((data['segment'] == 'NUVC')
                                   & (data['cenwave'] == cenwave))
                plot(ax, detector, cenwave, nuva_data=data[nuva_data], 
                     nuvb_data=data[nuvb_data], nuvc_data=data[nuvc_data])

        plt.tight_layout()
        filename = 'shifts_{}_{}.png'.format(detector, grating)
        remove_if_there(os.path.join(out_dir, filename))
        fig.savefig(os.path.join(out_dir, filename))
        plt.close(fig)
        os.chmod(os.path.join(out_dir, filename), 0o766)

def make_plots_per_cenwave(data, out_dir):
    """Plot shift vs time for different cenwaves.
    """
    sorted_index = np.argsort(data['date'])
    data = data[sorted_index]

    plt.rc('font', weight='bold')
    plt.rc('xtick.major', size=5, pad=7)
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15) 

    detectors = set(data['detector'])

    for detector in detectors:
        det_index = np.where(data['detector'] == detector)[0]
        cenwaves = set(data['cenwave'][det_index])

        for shift in ['x_shift', 'y_shift']:
            
            if shift == 'x_shift':
                shift_label = 'SHIFT1'
            else:
                shift_label = 'SHIFT2'

            for cenwave in cenwaves:
                if detector == 'FUV':
                    
                    f,(ax1,ax2) = plt.subplots(2,1,figsize=(30,15))
                    f.subplots_adjust(hspace=0.7)
                    # LP1
                    fuva_data_fp1 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 1) 
                                           & (data['life_adj'] ==1))[0]
                    
                    fuva_data_fp2 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 2)
                                           & (data['life_adj'] ==1))[0]
                    
                    fuva_data_fp3 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 3)
                                           & (data['life_adj'] ==1))[0]
                    
                    fuva_data_fp4 = np.where((data['segment'] == 'FUVA') 
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 4) 
                                           & (data['life_adj'] ==1))[0]
                    
                    fuvb_data_fp1 = np.where((data['segment'] == 'FUVB') 
                                           & (data['cenwave'] == cenwave)  
                                           & (data['fppos'] == 1) 
                                           & (data['life_adj'] ==1))[0]

                    fuvb_data_fp2 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave)  
                                           & (data['fppos'] == 2)
                                           & (data['life_adj'] ==1))[0]
                    
                    fuvb_data_fp3 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 3)
                                           & (data['life_adj'] ==1))[0]
                    
                    fuvb_data_fp4 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 4)
                                           & (data['life_adj'] ==1))[0]
                    
                    ax1.plot(data['date'][fuva_data_fp1], data[shift][fuva_data_fp1], 
                             'bo',label = 'FUVA LP1 | {} | FP1'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp2], data[shift][fuva_data_fp2], 
                             'bx',label = 'FUVA LP1 | {} | FP2'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp3], data[shift][fuva_data_fp3], 
                             'bv',label = 'FUVA LP1 | {} | FP3'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp4], data[shift][fuva_data_fp4], 
                             'b*',label = 'FUVA LP1 | {} | FP4'.format(cenwave))
                    
                    ax2.plot(data['date'][fuvb_data_fp1], data[shift][fuvb_data_fp1], 
                             'bo',label = 'FUVB LP1 | {} | FP1'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp2], data[shift][fuvb_data_fp2], 
                             'bx',label = 'FUVB LP1 | {} | FP2'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp3], data[shift][fuvb_data_fp3], 
                             'bv',label = 'FUVB LP1 | {} | FP3'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp4], data[shift][fuvb_data_fp4], 
                             'b*',label = 'FUVB LP1 | {} | FP4'.format(cenwave))
                    
                    # LP2
                    fuva_data_fp1 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 1) 
                                           & (data['life_adj'] ==2))[0]
                    
                    fuva_data_fp2 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 2)
                                           & (data['life_adj'] ==2))[0]
                    
                    fuva_data_fp3 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 3)
                                           & (data['life_adj'] ==2))[0]
                    
                    fuva_data_fp4 = np.where((data['segment'] == 'FUVA') 
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 4) 
                                           & (data['life_adj'] ==2))[0]
                    
                    fuvb_data_fp1 = np.where((data['segment'] == 'FUVB') 
                                           & (data['cenwave'] == cenwave)  
                                           & (data['fppos'] == 1) 
                                           & (data['life_adj'] ==2))[0]

                    fuvb_data_fp2 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave)  
                                           & (data['fppos'] == 2)
                                           & (data['life_adj'] ==2))[0]
                    
                    fuvb_data_fp3 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 3)
                                           & (data['life_adj'] ==2))[0]
                    
                    fuvb_data_fp4 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 4)
                                           & (data['life_adj'] ==2))[0]
                    
                    ax1.plot(data['date'][fuva_data_fp1], data[shift][fuva_data_fp1], 
                             'ro',label = 'FUVA LP2 | {} | FP1'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp2], data[shift][fuva_data_fp2], 
                             'rx',label = 'FUVA LP2 | {} | FP2'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp3], data[shift][fuva_data_fp3], 
                             'rv',label = 'FUVA LP2 | {} | FP3'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp4], data[shift][fuva_data_fp4], 
                             'r*',label = 'FUVA LP2 | {} | FP4'.format(cenwave))
                    
                    ax2.plot(data['date'][fuvb_data_fp1], data[shift][fuvb_data_fp1], 
                             'ro',label = 'FUVB LP2 | {} | FP1'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp2], data[shift][fuvb_data_fp2], 
                             'rx',label = 'FUVB LP2 | {} | FP2'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp3], data[shift][fuvb_data_fp3], 
                             'rv',label = 'FUVB LP2 | {} | FP3'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp4], data[shift][fuvb_data_fp4], 
                             'r*',label = 'FUVB LP2 | {} | FP4'.format(cenwave))
                    
                    # LP3
                    fuva_data_fp1 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 1) 
                                           & (data['life_adj'] ==3))[0]
                    
                    fuva_data_fp2 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 2)
                                           & (data['life_adj'] ==3))[0]
                    
                    fuva_data_fp3 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 3)
                                           & (data['life_adj'] ==3))[0]
                    
                    fuva_data_fp4 = np.where((data['segment'] == 'FUVA') 
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 4) 
                                           & (data['life_adj'] ==3))[0]
                    
                    fuvb_data_fp1 = np.where((data['segment'] == 'FUVB') 
                                           & (data['cenwave'] == cenwave)  
                                           & (data['fppos'] == 1) 
                                           & (data['life_adj'] ==3))[0]

                    fuvb_data_fp2 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave)  
                                           & (data['fppos'] == 2)
                                           & (data['life_adj'] ==3))[0]
                    
                    fuvb_data_fp3 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 3)
                                           & (data['life_adj'] ==3))[0]
                    
                    fuvb_data_fp4 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 4)
                                           & (data['life_adj'] ==3))[0]
                    
                    ax1.plot(data['date'][fuva_data_fp1], data[shift][fuva_data_fp1], 
                             'go',label = 'FUVA LP3 | {} | FP1'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp2], data[shift][fuva_data_fp2], 
                             'gx',label = 'FUVA LP3 | {} | FP2'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp3], data[shift][fuva_data_fp3], 
                             'gv',label = 'FUVA LP3 | {} | FP3'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp4], data[shift][fuva_data_fp4], 
                             'g*',label = 'FUVA LP3 | {} | FP4'.format(cenwave))
                    
                    ax2.plot(data['date'][fuvb_data_fp1], data[shift][fuvb_data_fp1], 
                             'go',label = 'FUVB LP3 | {} | FP1'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp2], data[shift][fuvb_data_fp2], 
                             'gx',label = 'FUVB LP3 | {} | FP2'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp3], data[shift][fuvb_data_fp3], 
                             'gv',label = 'FUVB LP3 | {} | FP3'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp4], data[shift][fuvb_data_fp4], 
                             'g*',label = 'FUVB LP3 | {} | FP4'.format(cenwave))
                    
                    #-- LP4
                    fuva_data_fp1 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 1) 
                                           & (data['life_adj'] ==4))[0]
                    
                    fuva_data_fp2 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 2)
                                           & (data['life_adj'] ==4))[0]
                    
                    fuva_data_fp3 = np.where((data['segment'] == 'FUVA')
                                           & (data['cenwave'] == cenwave)
                                           & (data['fppos'] == 3)
                                           & (data['life_adj'] ==4))[0]
                    
                    fuva_data_fp4 = np.where((data['segment'] == 'FUVA') 
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 4) 
                                           & (data['life_adj'] ==4))[0]
                    
                    fuvb_data_fp1 = np.where((data['segment'] == 'FUVB') 
                                           & (data['cenwave'] == cenwave)  
                                           & (data['fppos'] == 1) 
                                           & (data['life_adj'] ==4))[0]

                    fuvb_data_fp2 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave)  
                                           & (data['fppos'] == 2)
                                           & (data['life_adj'] ==4))[0]
                    
                    fuvb_data_fp3 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 3)
                                           & (data['life_adj'] ==4))[0]
                    
                    fuvb_data_fp4 = np.where((data['segment'] == 'FUVB')
                                           & (data['cenwave'] == cenwave) 
                                           & (data['fppos'] == 4)
                                           & (data['life_adj'] ==4))[0]
                    
                    ax1.plot(data['date'][fuva_data_fp1], data[shift][fuva_data_fp1], 
                             'ko',label = 'FUVA LP4 | {} | FP1'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp2], data[shift][fuva_data_fp2], 
                             'kx',label = 'FUVA LP4 | {} | FP2'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp3], data[shift][fuva_data_fp3], 
                             'kv',label = 'FUVA LP4 | {} | FP3'.format(cenwave))
                    ax1.plot(data['date'][fuva_data_fp4], data[shift][fuva_data_fp4], 
                             'k*',label = 'FUVA LP4 | {} | FP4'.format(cenwave))
                    
                    ax2.plot(data['date'][fuvb_data_fp1], data[shift][fuvb_data_fp1], 
                             'ko',label = 'FUVB LP4 | {} | FP1'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp2], data[shift][fuvb_data_fp2], 
                             'kx',label = 'FUVB LP4 | {} | FP2'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp3], data[shift][fuvb_data_fp3], 
                             'kv',label = 'FUVB LP4 | {} | FP3'.format(cenwave))
                    ax2.plot(data['date'][fuvb_data_fp4], data[shift][fuvb_data_fp4], 
                             'k*',label = 'FUVB LP4 | {} | FP4'.format(cenwave))
                    # End LPs

                    for lp_date in [56101,57062,58028]:
                        ax1.axvline(lp_date, ls='--', lw=2, c='k')
                        ax2.axvline(lp_date, ls='--', lw=2, c='k')

                    ax1.set_xlabel('DATE [MJD]', fontsize=20, 
                                   fontweight='bold')
                    ax1.set_ylabel('{} [Pixels]'.format(shift_label), 
                                   fontsize=20, fontweight='bold')
                    ax1.grid(linestyle='--')
                    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4),
                                ncol=4, fancybox=True, shadow=True, 
                                fontsize=15)
                    
                    ax2.set_xlabel('DATE [MJD]', fontsize=20, fontweight='bold')
                    ax2.set_ylabel('{} [Pixels]'.format(shift_label), 
                                   fontsize=20, fontweight='bold')
                    ax2.grid(linestyle='--')
                    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4),
                                ncol=4, fancybox=True, shadow=True, 
                                fontsize=15)

                    # f.suptitle('{}[A/B] {} | {}'.format(shift_label,detector,cenwave), fontsize=25, fontweight='bold')
                    filename = os.path.join(out_dir,
                                            'test','{}_{}_{}_BY_LIFE_POS.png'\
                                            .format(shift_label, 
                                                    detector, cenwave))
                    f.savefig(filename)
                    plt.close(f)

                elif detector == 'NUV':
                    f,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(30,15))
                    
                    #-- NUVA FP-POS 1-4
                    nuva_data_fp1 = np.where((data['segment'] == 'NUVA') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 1))[0]
                    
                    nuva_data_fp2 = np.where((data['segment'] == 'NUVA') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 2))[0]

                    nuva_data_fp3 = np.where((data['segment'] == 'NUVA') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 3))[0]

                    nuva_data_fp4 = np.where((data['segment'] == 'NUVA') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 4))[0]

                    #-- NUVB FP-POS 1-4
                    nuvb_data_fp1 = np.where((data['segment'] == 'NUVB') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 1))[0]
                    
                    nuvb_data_fp2 = np.where((data['segment'] == 'NUVB') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 2))[0]

                    nuvb_data_fp3 = np.where((data['segment'] == 'NUVB') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 3))[0]

                    nuvb_data_fp4 = np.where((data['segment'] == 'NUVB') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 4))[0]

                    #-- NUVC FP-POS 1-4
                    nuvc_data_fp1 = np.where((data['segment'] == 'NUVC') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 1))[0]
                    
                    nuvc_data_fp2 = np.where((data['segment'] == 'NUVC') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 2))[0]

                    nuvc_data_fp3 = np.where((data['segment'] == 'NUVC') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 3))[0]

                    nuvc_data_fp4 = np.where((data['segment'] == 'NUVC') &
                                            (data['cenwave'] == cenwave) &
                                            (data['fppos'] == 4))[0]

                    ax1.plot(data['date'][nuva_data_fp1], data[shift][nuva_data_fp1], 'bo',label = 'NUVA | {} | FP1'.format(cenwave))
                    ax1.plot(data['date'][nuva_data_fp2], data[shift][nuva_data_fp2], 'bx',label = 'NUVA | {} | FP2'.format(cenwave))
                    ax1.plot(data['date'][nuva_data_fp3], data[shift][nuva_data_fp3], 'bv',label = 'NUVA | {} | FP3'.format(cenwave))
                    ax1.plot(data['date'][nuva_data_fp4], data[shift][nuva_data_fp4], 'b*',label = 'NUVA | {} | FP4'.format(cenwave))
                    ax1.set_xlabel('DATE [MJD]', fontsize=20, fontweight='bold')
                    ax1.set_ylabel('{} [Pixels]'.format(shift_label), fontsize=20, fontweight='bold')
                    ax1.grid(linestyle='--')
                    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
                                ncol=4, fancybox=True, shadow=True, fontsize=15)

                    ax2.plot(data['date'][nuvb_data_fp1], data[shift][nuvb_data_fp1], 'ro',label = 'NUVB | {} | FP1'.format(cenwave))
                    ax2.plot(data['date'][nuvb_data_fp2], data[shift][nuvb_data_fp2], 'rx',label = 'NUVB | {} | FP2'.format(cenwave))
                    ax2.plot(data['date'][nuvb_data_fp3], data[shift][nuvb_data_fp3], 'rv',label = 'NUVB | {} | FP3'.format(cenwave))
                    ax2.plot(data['date'][nuvb_data_fp4], data[shift][nuvb_data_fp4], 'r*',label = 'NUVB | {} | FP4'.format(cenwave))
                    ax2.set_xlabel('DATE [MJD]', fontsize=20, fontweight='bold')
                    ax2.set_ylabel('{} [Pixels]'.format(shift_label), fontsize=20, fontweight='bold')
                    ax2.grid(linestyle='--')
                    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
                                ncol=4, fancybox=True, shadow=True, fontsize=15)

                    ax3.plot(data['date'][nuvc_data_fp1], data[shift][nuvc_data_fp1], 'go',label = 'NUVC | {} | FP1'.format(cenwave))
                    ax3.plot(data['date'][nuvc_data_fp2], data[shift][nuvc_data_fp2], 'gx',label = 'NUVC | {} | FP2'.format(cenwave))
                    ax3.plot(data['date'][nuvc_data_fp3], data[shift][nuvc_data_fp3], 'gv',label = 'NUVC | {} | FP3'.format(cenwave))
                    ax3.plot(data['date'][nuvc_data_fp4], data[shift][nuvc_data_fp4], 'g*',label = 'NUVC | {} | FP4'.format(cenwave))
                    ax3.set_xlabel('DATE [MJD]', fontsize=20, fontweight='bold')
                    ax3.set_ylabel('{} [Pixels]'.format(shift_label), fontsize=20, fontweight='bold')
                    ax3.grid(linestyle='--')
                    ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
                                ncol=4, fancybox=True, shadow=True, fontsize=15)

                    f.suptitle('{}[A/B/C] {} | {}'.format(shift_label,detector,cenwave), fontsize=25, fontweight='bold')
                    filename = os.path.join(out_dir,'{}_{}_{}.png'.format(shift_label, detector, cenwave))
                    f.savefig(filename)
                    plt.close(f)

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

def monitor():
    """Run the entire suite of monitoring
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    logger.info("STARTING MONITOR")

    settings = get_settings()

    webpage_dir = os.path.join(settings['webpage_location'], 'shifts')
    monitor_dir = os.path.join(settings['monitor_location'], 'Shifts')

    for place in [webpage_dir, monitor_dir]:
        if not os.path.exists(place):
            logger.debug("CREATING MONITOR LOCATION: {}".format(place))
            os.makedirs(place)

    # Create tables
    flash_data = make_shift_table(Lampflash)
    rawacq_data = make_shift_table(Rawacqs)
    
    # ascii.write(flash_data, os.path.join(monitor_dir,'monitor_db_table.csv'),
    #             format='csv', overwrite=True)

    # Make static plots.
    make_plots(flash_data, rawacq_data, monitor_dir)
    # make_plots_per_fppos(flash_data, monitor_dir)
    # make_plots_per_cenwave(flash_data, monitor_dir)
    make_interactive_plots(flash_data, rawacq_data, monitor_dir, 'FUV')
    make_interactive_plots(flash_data, rawacq_data, monitor_dir, 'NUV')
    # make_plots_2(flash_data, rawacq_data, monitor_dir)
    # fp_diff(flash_data)

    for item in glob.glob(os.path.join(monitor_dir, '*.p??')):
        remove_if_there(os.path.join(webpage_dir, os.path.basename(item)))
        shutil.copy(item, webpage_dir)

    logger.info("FINISH MONITOR")