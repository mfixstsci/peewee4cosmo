"""Script to monitor the COS FUV STIM pulses in TIME-TAG observations.

"""

from __future__ import absolute_import, print_function

__author__ = 'Justin Ely, Mees Fix'
__maintainer__ = 'Mees Fix'
__email__ = 'mfix@stsci.edu'
__status__ = 'Active'

import matplotlib.pyplot as plt
import datetime
import shutil
import os
import glob
import numpy as np
np.seterr(divide='ignore')
from astropy.table import Table
from astropy.time import Time
from astropy import table
from astropy.io import fits
from scipy.stats import linregress
import logging
logger = logging.getLogger(__name__)

from calcos import ccos

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

from ..database.models import get_database, get_settings
from ..database.models import Stims

from peewee import *

from ..utils import remove_if_there

import gc

from bokeh.io import output_file, show, save, gridplot
from bokeh.plotting import figure
from bokeh.models import Range1d, HoverTool, BoxSelectTool, ColumnDataSource, OpenURL, TapTool, Div, Button, CustomJS
from bokeh.layouts import column, row

from copy import deepcopy
#-------------------------------------------------------------------------------

def find_center(data):
    """ Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments.

    Parameters
    ----------
    data: np.array
        2d array of counts
    
    Returns
    -------
    x: int
        x position of the stim measurement
    y: int
        y position of the stim measurement
    """
    total = data.sum()
    if not total:
        return None, None

    X, Y = np.indices(data.shape)
    x = (X * data).sum() / total
    y = (Y * data).sum() / total
    #col = data[:, int(y)]
    #width_x = np.sqrt(abs((np.arange(col.size) - y) ** 2 * col).sum() / col.sum())
    #row = data[int(x), :]
    #width_y = np.sqrt(abs((np.arange(row.size) - x) ** 2 * row).sum() / row.sum())
    #height = data.max() 

    return y, x

#-------------------------------------------------------------------------------

def brf_positions(brftab, segment, position):
    """ Gets the search ranges for a stim pulse from
    the given baseline reference table.

    Parameters
    ----------
    brftab: str
        Path to brf reference file.
    segment: str
        FUVA of FUVB.
    position: str
        ul for upper left of lr for lower right.

    Returns
    -------
    xmin: int
        Minimum x vertex that defines the box the stims should fall in
    xmax: int
        Maximum x vertex that defines the box the stims should fall in
    ymin: int
        Mimimum y vertex that defines the box the stims should fall in
    ymax: int
        Maximum y vertex that defines the box the stims should fall in
    """
    brf = fits.getdata(brftab)

    if segment == 'FUVA':
        row = 0
    elif segment == 'FUVB':
        row = 1
    else:
        raise ValueError("Segment {} not understood.".format(segment))

    if position == 'ul':
        x_loc = 'SX1'
        y_loc = 'SY1'
    elif position == 'lr':
        x_loc = 'SX2'
        y_loc = 'SY2'

    xcenter = brf[row][x_loc]
    ycenter = brf[row][y_loc]
    xwidth = brf[row]['XWIDTH']
    ywidth = brf[row]['YWIDTH']

    xmin = int(max(xcenter - xwidth, 0))
    xmax = int(min(xcenter + xwidth, 16384))
    ymin = int(max(ycenter - ywidth, 0))
    ymax = int(min(ycenter + ywidth, 1024))

    return xmin, xmax, ymin, ymax

#-------------------------------------------------------------------------------

def find_stims(image, segment, stim, brf_file):
    """Locate stims for table insertion.

    Parameters
    ----------
    image: np.array
        A numpy array of counts.
    segment: str
        FUVA or FUVB
    stim: str
        ul for upper left or lr for lower right
    brf_file: str
        Path to brf reference file.
    """
    x1, x2, y1, y2 = brf_positions(brf_file, segment, stim)
    found_x, found_y = find_center(image[y1:y2, x1:x2])
    if (not found_x) and (not found_y):
        return -999, -999

    return found_x + x1, found_y + y1

#-------------------------------------------------------------------------------

def locate_stims(data_object, start=0, increment=None, brf_file=None):
    """Breaks a timetag exposure into different time steps and measures
    the x + y poisitons of the stim pulses to track them over the 
    length of an exposure.

    Parameters
    ----------
    data_object: peewee row object
        Row from the cosmo database with path and filename attributes.
    start: int
        Time step start
    increment: float
        Time step increment (dt)
    brf_file: str
        path and filename for brf reference file.

    Yields
    ------
    info: dict
        Dictionary with information that will make a row in the DB.
    """
    full_path = os.path.join(data_object.path, data_object.filename)

    #-- change this to pull brf file from the header if not specified    
    if not brf_file:
        brf_file = os.path.join(os.environ['lref'], 's7g1700el_brf.fits')

    DAYS_PER_SECOND = 1. / 60. / 60. / 24.

    file_path, file_name = os.path.split(full_path)

    with fits.open(full_path) as hdu:
        exptime = hdu[1].header['exptime']
        expstart = hdu[1].header['expstart']
        segment = hdu[0].header['segment']

        stim_info = {'filename': file_name,
                     'rootname': hdu[0].header['rootname'],
                     'proposid': hdu[0].header['proposid'],
                     'segment': segment}

        try:
            hdu[1].data
        except:
            yield stim_info
            raise StopIteration

        if not len(hdu[1].data):
            yield stim_info
            raise StopIteration

        #-- If increment is not supplied, use the rates supplied by the detector
        if not increment:
            if exptime < 10:
                increment = .03
            elif exptime < 100:
                increment = 1
            else:
                increment = 30

            increment *= 4

        stop = start + increment

        #-- Iterate from start to stop, excluding final bin if smaller than increment
        start_times = np.arange(start, exptime-increment, increment)
        if not len(start_times):
            yield stim_info

        for sub_start in start_times:
            events = hdu['events'].data

            #-- No COS observation has data below ~923
            data_index = np.where((hdu[1].data['time'] >= sub_start) &
                                  (hdu[1].data['time'] <= sub_start+increment))[0]

            events = events[data_index]

            #-- Call for this is x_values, y_values, image to bin to, offset in x
            # ccos.binevents(x, y, array, x_offset, dq, sdqflags, epsilon)
            im = np.zeros((1024, 16384)).astype(np.float32)
            ccos.binevents(events['RAWX'].astype(np.float32),
                           events['RAWY'].astype(np.float32),
                           im,
                           0,
                           events['dq'],
                           0)

            ABS_TIME = expstart + sub_start * DAYS_PER_SECOND

            found_ul_x, found_ul_y = find_stims(im, segment, 'ul', brf_file)
            found_lr_x, found_lr_y = find_stims(im, segment, 'lr', brf_file)

            stim_info['time'] = round(sub_start, 5)
            stim_info['abs_time'] = round(ABS_TIME, 5)
            stim_info['stim1_x'] = round(found_ul_x, 3)
            stim_info['stim1_y'] = round(found_ul_y, 3)
            stim_info['stim2_x'] = round(found_lr_x, 3)
            stim_info['stim2_y'] = round(found_lr_y, 3)
            stim_info['counts'] = round(im.sum(), 7)

            yield deepcopy(stim_info)
        
        #-- delete hdu at the end
        #-- collect the garbage.
        del hdu
        gc.collect()

#-------------------------------------------------------------------------------

def make_position_panel(segment, stim):
    
    """Make a bokeh plotting panel for the x + y 
    pixel positions.

    Parameters
    ----------
    segment: str
        FUVA or FUVB
    stim: str
        upper left or lower right.

    Returns
    -------
    p: bokeh panel
        Panel containing scatter plot for stim/segment combo.
    """
    
    #-- Open brftab and get some data.
    brf_file = os.path.join(os.environ['lref'], 's7g1700el_brf.fits')
    brf = fits.getdata(brf_file, 1)

    #-- Build "box" that stims should be contained in.
    if segment == 'FUVA' and stim == 'upper':
        xcenter = brf[0]['SX1']
        ycenter = brf[0]['SY1']
        xwidth = brf[0]['XWIDTH']
        ywidth = brf[0]['YWIDTH']

        x_range_min, x_range_max = 250, 530
        y_range_min, y_range_max = 940, 1030
    elif segment == 'FUVA' and stim == 'lower':
        xcenter = brf[0]['SX2']
        ycenter = brf[0]['SY2']
        xwidth = brf[0]['XWIDTH']
        ywidth = brf[0]['YWIDTH']
        
        x_range_min, x_range_max = 15860, 16140
        y_range_min, y_range_max = -15, 80
    elif segment == 'FUVB' and stim == 'upper':
        xcenter = brf[1]['SX1']
        ycenter = brf[1]['SY1']
        xwidth = brf[1]['XWIDTH']
        ywidth = brf[1]['YWIDTH']

        x_range_min, x_range_max = 320, 460
        y_range_min, y_range_max = 875, 1100
    elif segment == 'FUVB' and stim == 'lower':
        xcenter = brf[1]['SX2']
        ycenter = brf[1]['SY2']
        xwidth = brf[1]['XWIDTH']
        ywidth = brf[1]['YWIDTH']

        x_range_min, x_range_max = 15920, 16080
        y_range_min, y_range_max = -75, 150
    else:
        logger.info('WE DONT SUPPORT CONFIG {} {}'.formar(segment, stim))
    
    #-- Make box where stims should ideally be...
    xs = [xcenter - xwidth,
          xcenter + xwidth,
          xcenter + xwidth,
          xcenter - xwidth,
          xcenter - xwidth]
    ys = [ycenter - ywidth,
          ycenter - ywidth,
          ycenter + ywidth,
          ycenter + ywidth,
          ycenter - ywidth]

    #-- Connect to DB
    database = get_database()
    database.connect()

    #-- Set height and width.
    plt_wth = 600
    plt_hgt = 500

    if segment == 'FUVA':
        color = 'blue'
    else:
        color = 'red'

    TOOLS ='box_zoom,box_select,crosshair,pan,reset,tap,save'

    #data = Stims.select().where(Stims.segment==segment)
    data = []
    for i, row in enumerate(Stims.select().where(Stims.segment==segment).dicts()):
        data.append(row.values())
        if not i:
            #-- get keys here because if you use ._meta.fields.keys() 
            #-- they will be out of order.
            keys = row.keys()    
    
    database.close()

    data = Table(rows=data, names=keys)
    data = table.unique(data,['rootname',
                              'stim1_x',
                              'stim2_x',
                              'stim1_y',
                              'stim2_y'])
    if stim == 'upper':
        #-- Build ColumnDataSource object
        source = ColumnDataSource(data=dict(
                        date=data['time'],
                        rootname=data['rootname'],
                        stim1_x=data['stim1_x'],
                        stim1_y=data['stim1_y'],
                        segment=data['segment'],
                        abs_time=data['abs_time'],
                        counts=data['counts'],
                        proposid=data['proposid'],
                        ))
        
        stim_x, stim_y = 'stim1_x', 'stim1_y'
        location = 'Upper Left'
        
    else:
        source = ColumnDataSource(data=dict(
                        date=data['time'],
                        rootname=data['rootname'],
                        stim2_x=data['stim2_x'],
                        stim2_y=data['stim2_y'],
                        segment=data['segment'],
                        abs_time=data['abs_time'],
                        counts=data['counts'],
                        proposid=data['proposid'],
                        ))
        
        stim_x, stim_y = 'stim2_x', 'stim2_y'
        location = 'Lower Right'
        
    #-- Build Hovertool Tips.
    hover = HoverTool(
                tooltips=[
                        ("(x, y):","($x, $y)"),
                        ("Date", "@abs_time"),
                        ("Rootname","@rootname"),
                        ("Segment","@segment"),
                        ("Counts","@counts"),
                        ("Proposal ID","@proposid")
                        ],
                names=['data_scatter']   
                     )

    p = figure(width=plt_wth, height=plt_hgt, x_range=Range1d(x_range_min, x_range_max), y_range=Range1d(y_range_min, y_range_max), title='{} {} Stim'.format(segment, location),tools=[TOOLS,hover])
    p.title.text_font_size = '15pt'
    p.circle(stim_x, stim_y, legend='Stim Pulse Location',size=4, name='data_scatter', color="black", source=source, alpha=0.5)
    p.line(xs, ys, legend='Box Boundary', color=color, line_width=4, line_dash='dashed')
    p.xaxis.axis_label = 'RAWX (Pixels)'
    p.yaxis.axis_label = 'RAWY (Pixels)'
    p.xaxis.axis_label_text_font_size = "13pt"
    p.yaxis.axis_label_text_font_size = "13pt"

    #-- Provide URL and taptool and callback info.
    url = "http://archive.stsci.edu/proposal_search.php?id=@proposid&mission=hst"
    taptool = p.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

    return p
#-------------------------------------------------------------------------------
def make_time_plot(segment):
    """Make the plot of the stim positions in 1d as a function of time for x + y
    for each segment.

    Parameters
    ----------
    segment: str
        FUVA or FUVB

    Returns
    -------
    None
    """
    settings = get_settings()
    monitor_dir = settings['monitor_location']

    outname = os.path.join(monitor_dir, 'stim_time_{}.html'.format(segment))
    remove_if_there(outname)
    output_file(outname)

    #-- Connect to DB
    database = get_database()
    database.connect()

    #-- Set height and width.
    plt_wth = 600
    plt_hgt = 500

    #data = Stims.select().where(Stims.segment==segment)
    col_names = ['stim1_x', 'stim1_y', 'stim2_x', 'stim2_y']
    titles = ['Upper Left, X', 'Upper Left, Y', 'Lower Right, X', 'Lower Right, Y']
    
    ylims = {'FUVA':{'stim1_x': [390, 425],
                     'stim1_y': [950, 990],
                     'stim2_x': [16000, 16120],
                     'stim2_y': [-10, 40]},
             'FUVB':{'stim1_x': [370, 430],
                     'stim1_y': [920, 1000],
                     'stim2_x': [15960, 16070],
                     'stim2_y': [-10, 40]}}


    panels = []
    for stim_loc, stim in enumerate(col_names):
        query = list(Stims.raw("""SELECT * 
                                    FROM stims 
                                   WHERE segment = %s""", segment).dicts())

        time = [row['abs_time'] for row in query]
        stims = [row[stim] for row in query]

        p = figure(width=plt_wth, height=plt_hgt, x_range=Range1d(54900, max(time)+100), y_range=Range1d(ylims[segment][stim][0], ylims[segment][stim][1]), title='{} {} Stim'.format(segment, titles[stim_loc]))
        p.title.text_font_size = '15pt'
        p.circle(time, stims, legend='Stim Pulse Location',size=4, name='data_scatter', color="black", alpha=0.5)
        p.xaxis.axis_label = 'Time (MJD)'
        p.yaxis.axis_label = titles[stim_loc] + ' (Pixels)'
        p.xaxis.axis_label_text_font_size = "13pt"
        p.yaxis.axis_label_text_font_size = "13pt"

        panels.append(p)
    
    database.close()
    
    full_plot = gridplot([[panels[0], panels[1]], [panels[2], panels[3]]])
    save(full_plot, filename=outname)
#-------------------------------------------------------------------------------
def make_stretch_panel(segment, left_stim, right_stim, reverse_stims=False):
    """Make plotting panel of stim 'stretch' in x or y.

    Parameters
    ----------
    segment: str
        FUVA or FUVB
    left_stim: str
        stim of lower coordinate value
    right_stim: str
        stim of higher coordinate value
    reverse_stims: bool
        reverse the order the stims are subtracted from each other
    
    Returns
    -------
    p: bokeh panel
        Panel containing scatter of stretch vs. time
    """
    #-- Connect to DB
    database = get_database()
    database.connect()

    #-- Set height and width.
    plt_wth = 600
    plt_hgt = 500

    #data = Stims.select().where(Stims.segment==segment)
    
    query = list(Stims.select().where(Stims.segment==segment).dicts())

    database.close()
    
    time = [row['abs_time'] for row in query]

    if reverse_stims:
        stretch = [row[left_stim] - row[right_stim] for row in query]
        ylims = {'FUVA': [945, 980],
                 'FUVB': [930, 980]}
        label_str = '(Stim 1 y) - (Stim 2 y)'
    else:
        stretch = [row[right_stim] - row[left_stim] for row in query]
        ylims = {'FUVA': [15600,15740],
                 'FUVB': [15560, 15700]}
        label_str = '(Stim 2 x) - (Stim 1 x)'

    p = figure(width=plt_wth, height=plt_hgt, x_range=Range1d(54900, max(time)+100), y_range=Range1d(ylims[segment][0], ylims[segment][1]), title='{} {}'.format(segment, label_str))
    p.title.text_font_size = '15pt'
    p.circle(time, stretch, legend='Stretch (Pixels)',size=4, name='data_scatter', color="black", alpha=0.5)
    p.xaxis.axis_label = 'Time (MJD)'
    p.yaxis.axis_label = 'Stretch'
    p.xaxis.axis_label_text_font_size = "13pt"
    p.yaxis.axis_label_text_font_size = "13pt"

    return p
#-------------------------------------------------------------------------------

def interactive_plotting(path=None, filename=None, position=False, time=False, stretch=False):
    """Make interactive Bokeh Figures

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    if position:
        outname = os.path.join(path, filename)
        remove_if_there(outname)
        output_file(outname) 
        
        #-- FUV Panels x vs. y
        FUVA_upper = make_position_panel('FUVA', 'upper')
        FUVA_lower = make_position_panel('FUVA', 'lower')
        FUVB_upper = make_position_panel('FUVB', 'upper')
        FUVB_lower = make_position_panel('FUVB', 'lower')
        
        p = gridplot([[FUVA_upper, FUVA_lower], [FUVB_upper, FUVB_lower]])

        save(p, filename=outname)
    
    if time:
        #-- FUV stim x and y pos vs time.
        make_time_plot('FUVA')
        make_time_plot('FUVB')  

    if stretch:
        outname = os.path.join(path, filename)
        remove_if_there(outname)
        output_file(outname) 
        
        #-- FUV Panels x vs. y
        FUVA_x_stretch = make_stretch_panel('FUVA', 'stim1_x', 'stim2_x')
        FUVA_y_stretch = make_stretch_panel('FUVA', 'stim1_y', 'stim2_y', reverse_stims=True)
        FUVB_x_stretch = make_stretch_panel('FUVB', 'stim1_x', 'stim2_x')
        FUVB_y_stretch = make_stretch_panel('FUVB', 'stim1_y', 'stim2_y', reverse_stims=True)
        
        p = gridplot([[FUVA_x_stretch, FUVA_y_stretch], [FUVB_x_stretch, FUVB_y_stretch]])

        save(p, filename=outname)
#-------------------------------------------------------------------------------

def stim_monitor():
    """Main function to monitor the stim pulses in COS observations

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    logger.info("STARTING MONITOR")

    settings = get_settings()

    webpage_dir = os.path.join(settings['webpage_location'], 'stim')
    monitor_dir = os.path.join(settings['monitor_location'], 'Stims')
    
    logger.info("MAKING INTERACTIVE PLOTS")
    interactive_plotting(path=monitor_dir, filename='stim_location.html', position=True)
    logger.info("POSITION DONE")
    interactive_plotting(time=True)
    logger.info("TIME DONE")
    interactive_plotting(path=monitor_dir, filename='stim_stretch.html', stretch=True)
    logger.info("STRETCH DONE")
    # move_to_web(monitor_dir, webpage_dir)
    logger.info("FINISH MONITOR")