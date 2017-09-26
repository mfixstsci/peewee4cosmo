"""Script to monitor the COS FUV STIM pulses in TIME-TAG observations.

"""

from __future__ import absolute_import, print_function

__author__ = 'Justin Ely, Mees Fix'
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
    moments
    """
    total = data.sum()
    if not total:
        return None, None

    X, Y = np.indices(data.shape)
    x = (X * data).sum() // total
    y = (Y * data).sum() // total
    #col = data[:, int(y)]
    #width_x = np.sqrt(abs((np.arange(col.size) - y) ** 2 * col).sum() / col.sum())
    #row = data[int(x), :]
    #width_y = np.sqrt(abs((np.arange(row.size) - x) ** 2 * row).sum() / row.sum())
    #height = data.max()

    return y, x

#-------------------------------------------------------------------------------

def brf_positions(brftab, segment, position):
    """ Gets the search ranges for a stimpulse from
    the given baseline reference table
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
    x1, x2, y1, y2 = brf_positions(brf_file, segment, stim)
    found_x, found_y = find_center(image[y1:y2, x1:x2])
    if (not found_x) and (not found_y):
        return -999, -999

    return found_x + x1, found_y + y1

#-------------------------------------------------------------------------------

def locate_stims(data_object, start=0, increment=None, brf_file=None):

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
#-----------------------------------------------------

def find_missing():
    """Return list of datasets that had a missing stim in
    any time period.
    """

    SETTINGS = open_settings()
    Session, engine = load_connection(SETTINGS['connection_string'])

    data = engine.execute("""SELECT headers.rootname,stims.abs_time
                                    FROM stims
                                    JOIN headers ON stims.rootname = headers.rootname
                                    WHERE stims.stim1_x = -999
                                        OR stims.stim1_y = -999
                                        OR stims.stim2_x = -999
                                        OR stims.stim2_y = -999
                                    ORDER BY stims.abs_time""")

    missed_data = [(item[0].strip(), float(item[1])) for item in data]
    all_obs = [item[0] for item in missed_data]
    all_mjd = [item[1] for item in missed_data]

    obs_uniq = list(set(all_obs))

    date_uniq = []
    for obs in obs_uniq:
        min_mjd = all_mjd[all_obs.index(obs)]
        date_uniq.append(min_mjd)

    date_uniq = Time(date_uniq, format='mjd', scale='utc')
    return obs_uniq, date_uniq

#-------------------------------------------------------------------------------

def check_individual(out_dir, connection_string):
    """Run tests on each individual datasets.

    Currently checks if any coordinate deviates from the mean
    over the exposure and produces a plot if so.

    """

    Session, engine = load_connection(connection_string)

    query = """SELECT headers.rootname,
                      headers.segment,
                      headers.proposid,
                      headers.targname,
                      STD(stims.stim1_x) as stim1_xstd,
                      STD(stims.stim1_y) as stim1_ystd,
                      STD(stims.stim2_x) as stim2_xstd,
                      STD(stims.stim2_y) as stim2_ystd
                      FROM stims
                      JOIN headers on stims.rootname = headers.rootname
                      GROUP BY stims.file_id,
                               headers.rootname,
                               headers.segment,
                               headers.proposid,
                               headers.targname
                      HAVING stim1_xstd > 2 OR
                             stim1_ystd > 2 OR
                             stim2_xstd > 2 OR
                             stim2_ystd > 2;"""

    data = []
    for i, row in enumerate(engine.execute(query)):
        if not i:
            keys = row.keys()
        data.append(row.values())

    t = Table(rows=data, names=keys)
    t.write(os.path.join(out_dir, "STIM_problem_rootnames.txt"), delimiter='|', format='ascii.fixed_width')

#-------------------------------------------------------------------------------
def make_panel(segment, stim, time=False):
    

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
    plt_wth = 500
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

def interactive_plotting(path, filename):
    """Make interactive Bokeh Figures

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    outname = os.path.join(path, filename)
    remove_if_there(outname)
    output_file(outname) 
    
    #-- FUVA Panels x vs. y
    FUVA_upper = make_panel('FUVA', 'upper')
    FUVA_lower = make_panel('FUVA', 'lower')
    FUVB_upper = make_panel('FUVB', 'upper')
    FUVB_lower = make_panel('FUVB', 'lower')
    
    p = gridplot([[FUVA_upper, FUVA_lower], [FUVB_upper, FUVB_lower]])

    save(p, filename=outname)
#-------------------------------------------------------------------------------

def make_plots(out_dir, connection_string):
    """Make the overall STIM monitor plots.
    They will all be output to out_dir.
    """

    plt.ioff()

    brf_file = os.path.join(os.environ['lref'], 's7g1700el_brf.fits')

    brf = fits.getdata(brf_file, 1)

    database = get_database()
    database.connect()

    plt.figure(1, figsize=(18, 12))
    plt.grid(True)

    #-- stim1_x, stim1_y FUVA
    # data = engine.execute("""SELECT stim1_x, stim1_y
    #                                 FROM stims
    #                                 JOIN headers on stims.rootname = headers.rootname
    #                                 WHERE headers.segment = 'FUVA' AND
    #                                     stims.stim1_x != -999 AND
    #                                     stims.stim1_y != -999 AND
    #                                     stims.stim2_x != -999 AND
    #                                     stims.stim2_y != -999;""")

    data = Stims.select().where(Stims.segment=='FUVA')
    
    data = [line for line in data]
    plt.subplot(2, 2, 1)

    #-- FUVA plots
    x_range_min, x_range_max = 250, 530
    y_range_min, y_range_max = 940, 1030

    x = [line.stim1_x for line in data]
    y = [line.stim1_y for line in data]
    plt.plot(x, y, 'b.', alpha=.7)
    
    plt.xlim([x_range_min, x_range_max])
    plt.ylim([y_range_min, y_range_max])

    plt.xlabel('x')
    plt.ylabel('y')
    
    plt.title('Segment A: Stim A (Upper Left)')
    xcenter = brf[0]['SX1']
    ycenter = brf[0]['SY1']
    xwidth = brf[0]['XWIDTH']
    ywidth = brf[0]['YWIDTH']
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
    plt.plot(xs, ys, color='r', linestyle='--', label='Search Box')
    plt.legend(shadow=True, numpoints=1)
    plt.xlabel('RAWX')
    plt.ylabel('RAWY')
    #plt.set_xlims(xcenter - 2*xwidth, xcenter + 2*xwidth)
    #plt.set_ylims(ycenter - 2*ywidth, ycenter - 2*ywidth)

    #-- stim2_x, stim2_y FUVA
    # data = engine.execute("""SELECT stim2_x, stim2_y
    #                                 FROM stims
    #                                 JOIN headers on stims.rootname = headers.rootname
    #                                 WHERE headers.segment = 'FUVA' AND
    #                                     stims.stim1_x != -999 AND
    #                                     stims.stim1_y != -999 AND
    #                                     stims.stim2_x != -999 AND
    #                                     stims.stim2_y != -999;""")

    #-- FUVA Lower
    data = Stims.select().where(Stims.segment=='FUVA')
    data = [line for line in data]
   
    x_range_min, x_range_max = 15860, 16140
    y_range_min, y_range_max = -15, 80
   
    plt.subplot(2, 2, 2)
    x = [line.stim2_x for line in data]
    y = [line.stim2_y for line in data]
    plt.plot(x, y, 'r.', alpha=.7)

    plt.xlim([x_range_min, x_range_max])
    plt.ylim([x_range_min, x_range_max])

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Segment A: Stim B (Lower Right)')
    xcenter = brf[0]['SX2']
    ycenter = brf[0]['SY2']
    xwidth = brf[0]['XWIDTH']
    ywidth = brf[0]['YWIDTH']
    xs = [xcenter - xwidth,
          xcenter + xwidth,
          xcenter + xwidth,
          xcenter - xwidth,
          xcenter - xwidth]
    ys = [ycenter - ywidth,
          ycenter  -ywidth,
          ycenter + ywidth,
          ycenter + ywidth,
          ycenter - ywidth]
    plt.plot(xs, ys, color='r', linestyle='--', label='Search Box')
    plt.legend(shadow=True, numpoints=1)
    plt.xlabel('RAWX')
    plt.ylabel('RAWY')

    #-- stim1_x stim1_y FUVB
    # data = engine.execute("""SELECT stim1_x, stim1_y
    #                                 FROM stims
    #                                 JOIN headers on stims.rootname = headers.rootname
    #                                 WHERE headers.segment = 'FUVB' AND
    #                                     stims.stim1_x != -999 AND
    #                                     stims.stim1_y != -999 AND
    #                                     stims.stim2_x != -999 AND
    #                                     stims.stim2_y != -999;""")

    #-- FUVB Upper
    data = Stims.select().where(Stims.segment=='FUVB')
    data = [line for line in data]

    x_range_min, x_range_max = 320, 460
    y_range_min, y_range_max = 875, 1100
    
    plt.subplot(2, 2, 3)
    x = [line.stim1_x for line in data]
    y = [line.stim1_y for line in data]
    plt.plot(x, y, 'b.', alpha=.7)
    
    plt.xlim([x_range_min, x_range_max])
    plt.ylim([y_range_min, y_range_max])
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Segment B: Stim A (Upper Left)')
    xcenter = brf[1]['SX1']
    ycenter = brf[1]['SY1']
    xwidth = brf[1]['XWIDTH']
    ywidth = brf[1]['YWIDTH']
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
    plt.plot(xs, ys, color='r', linestyle='--', label='Search Box')
    plt.legend(shadow=True, numpoints=1)
    plt.xlabel('RAWX')
    plt.ylabel('RAWY')

    #-- stim2_x stim2_y FUVB
    # data = engine.execute("""SELECT stim2_x, stim2_y
    #                                 FROM stims
    #                                 JOIN headers on stims.rootname = headers.rootname
    #                                 WHERE headers.segment = 'FUVB' AND
    #                                     stims.stim1_x != -999 AND
    #                                     stims.stim1_y != -999 AND
    #                                     stims.stim2_x != -999 AND
    #                                     stims.stim2_y != -999;""")
    
    #-- FUVB Lower
    data = Stims.select().where(Stims.segment=='FUVB')
    data = [line for line in data]

    plt.subplot(2, 2, 4)
    
    x_range_min, x_range_max = 15920, 16080
    y_range_min, y_range_max = -75, 150
    
    x = [line.stim2_x for line in data]
    y = [line.stim2_y for line in data]
    plt.plot(x, y, 'r.', alpha=.7)

    plt.xlim([x_range_min, x_range_max])
    plt.ylim([y_range_min, y_range_max])

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Segment B: Stim B (Lower Right)')
    xcenter = brf[1]['SX2']
    ycenter = brf[1]['SY2']
    xwidth = brf[1]['XWIDTH']
    ywidth = brf[1]['YWIDTH']
    xs = [xcenter - xwidth,
          xcenter +xwidth,
          xcenter + xwidth,
          xcenter - xwidth,
          xcenter - xwidth]
    ys = [ycenter - ywidth,
          ycenter - ywidth,
          ycenter + ywidth,
          ycenter + ywidth,
          ycenter - ywidth]
    plt.plot(xs, ys, color='r', linestyle='--', label='Search Box')
    plt.legend(shadow=True, numpoints=1)
    plt.xlabel('RAWX')
    plt.ylabel('RAWY')

    plt.draw()
    remove_if_there(os.path.join(out_dir, 'STIM_locations_peewee.png'))
    plt.savefig(os.path.join(out_dir, 'STIM_locations_peewee.png'))
    plt.close(1)
    os.chmod(os.path.join(out_dir, 'STIM_locations_peewee.png'),0o766)


    # for segment in ['FUVA', 'FUVB']:
    #     fig = plt.figure(2, figsize=(18, 12))
    #     fig.suptitle('%s coordinate locations with time' % (segment))

    #     col_names = ['stim1_x', 'stim1_y', 'stim2_x', 'stim2_y']
    #     titles = ['Upper Left, X', 'Upper Left, Y', 'Lower Right, X', 'Lower Right, Y']

    #     for i, (column, title) in enumerate(zip(col_names, titles)):
    #         ax = fig.add_subplot(2, 2, i + 1)
    #         ax.set_title(title)
    #         ax.set_xlabel('MJD')
    #         ax.set_ylabel('Coordinate')

    #         #-- Data for stim vs time.
    #         # query = """SELECT stims.abs_time, stims.{}
    #         #                   FROM stims
    #         #                   JOIN headers ON stims.rootname = headers.rootname
    #         #                   WHERE headers.segment = '{}' AND
    #         #                       stims.stim1_x != -999 AND
    #         #                       stims.stim1_y != -999 AND
    #         #                       stims.stim2_x != -999 AND
    #         #                       stims.stim2_y != -999;""".format(column, segment)

    #         query = Stims.raw('SELECT abs_time, %s FROM stims WHERE segment = %s', (column,segment))
            
    #         data = [line for line in database.execute_sql(query)]
    #         times = [line[0] for line in data]
    #         coords = [line[1] for line in data]
    #         ax.plot(times, coords, 'o')

    #     remove_if_there(os.path.join(out_dir, 'STIM_locations_vs_time_%s_peewee.png' %
    #                                                                 (segment)))
    #     fig.savefig(os.path.join(out_dir, 'STIM_locations_vs_time_%s_peewee.png' %
    #                                                                 (segment)))
    #     plt.close(fig)
    #     os.chmod(os.path.join(out_dir, 'STIM_locations_vs_time_%s_peewee.png' %
    #                                                                 (segment)),0o766)

    database.close()
    # ------------------------#
    # strech and midpoint     #
    # ------------------------#
    # for segment in ['FUVA', 'FUVB']:
    #     fig = plt.figure(figsize=(18, 12))
    #     fig.suptitle("Strech and Midpoint vs time")

    #     ax1 = fig.add_subplot(2, 2, 1)
    #     query = """SELECT stims.abs_time, stims.stim2_x - stims.stim1_x as stretch
    #                       FROM stims
    #                       JOIN headers ON stims.rootname = headers.rootname
    #                       WHERE headers.segment = '{}' AND
    #                           stims.stim1_x != -999 AND
    #                           stims.stim1_y != -999 AND
    #                           stims.stim2_x != -999 AND
    #                           stims.stim2_y != -999;""".format(segment)
    #     data = [line for line in engine.execute(query)]
    #     stretch = [line.stretch for line in data]
    #     times = [line.abs_time for line in data]

    #     ax1.plot(times, stretch, 'o')
    #     ax1.set_xlabel('MJD')
    #     ax1.set_ylabel('Stretch X')

    #     ax2 = fig.add_subplot(2, 2, 2)
    #     query = """SELECT stims.abs_time, .5*(stims.stim2_x + stims.stim1_x) as midpoint
    #                       FROM stims
    #                       JOIN headers ON stims.rootname = headers.rootname
    #                       WHERE headers.segment = '{}' AND
    #                           stims.stim1_x != -999 AND
    #                           stims.stim1_y != -999 AND
    #                           stims.stim2_x != -999 AND
    #                           stims.stim2_y != -999;""".format(segment)
    #     data = [line for line in engine.execute(query)]
    #     midpoint = [line.midpoint for line in data]
    #     times = [line.abs_time for line in data]

    #     ax2.plot(times, midpoint, 'o')
    #     ax2.set_xlabel('MJD')
    #     ax2.set_ylabel('Midpoint X')

    #     ax3 = fig.add_subplot(2, 2, 3)
    #     query = """SELECT stims.abs_time, stims.stim2_y - stims.stim1_y as stretch
    #                       FROM stims
    #                       JOIN headers ON stims.rootname = headers.rootname
    #                       WHERE headers.segment = '{}' AND
    #                           stims.stim1_x != -999 AND
    #                           stims.stim1_y != -999 AND
    #                           stims.stim2_x != -999 AND
    #                           stims.stim2_y != -999;""".format(segment)
    #     data = [line for line in engine.execute(query)]
    #     stretch = [line.stretch for line in data]
    #     times = [line.abs_time for line in data]
    #     ax3.plot(times, stretch, 'o')
    #     ax3.set_xlabel('MJD')
    #     ax3.set_ylabel('Stretch Y')

    #     ax4 = fig.add_subplot(2, 2, 4)
    #     query = """SELECT stims.abs_time, .5*(stims.stim2_y + stims.stim1_y) as midpoint
    #                       FROM stims
    #                       JOIN headers ON stims.rootname = headers.rootname
    #                       WHERE headers.segment = '{}' AND
    #                           stims.stim1_x != -999 AND
    #                           stims.stim1_y != -999 AND
    #                           stims.stim2_x != -999 AND
    #                           stims.stim2_y != -999;""".format(segment)
    #     ax4.plot(times, midpoint, 'o')
    #     ax4.set_xlabel('MJD')
    #     ax4.set_ylabel('Midpoint Y')
    #     remove_if_there(os.path.join(out_dir, 'STIM_stretch_vs_time_%s.png' %
    #                  (segment)))
    #     fig.savefig(
    #         os.path.join(out_dir, 'STIM_stretch_vs_time_%s.png' %
    #                      (segment)))
    #     plt.close(fig)
    #     os.chmod(os.path.join(out_dir, 'STIM_stretch_vs_time_%s.png' %
    #                  (segment)),0o766)

    # fig = plt.figure(1, figsize=(18, 12))
    # ax = fig.add_subplot(2, 2, 1)
    # ax.grid(True)

    # data = engine.execute("""SELECT stim1_x, stim2_x
    #                                 FROM stims
    #                                 JOIN headers on stims.rootname = headers.rootname
    #                                 WHERE headers.segment = 'FUVA' AND
    #                                     stims.stim1_x != -999 AND
    #                                     stims.stim1_y != -999 AND
    #                                     stims.stim2_x != -999 AND
    #                                     stims.stim2_y != -999;""")
    # data = [line for line in data]

    # x1 = [float(line.stim1_x) for line in data]
    # x2 = [float(line.stim2_x) for line in data]

    # im, nothin1, nothin2 = np.histogram2d(x2, x1, bins=200)  ##reverse coords
    # im = np.log(im)
    # ax.imshow(im, aspect='auto', interpolation='none')
    # ax.set_xlabel('x1')
    # ax.set_ylabel('x2')
    # ax.set_title('Segment A: X vs X')



    # ax = fig.add_subplot(2, 2, 2)
    # ax.grid(True)

    # data = engine.execute("""SELECT stim1_y, stim2_y
    #                                 FROM stims
    #                                 JOIN headers on stims.rootname = headers.rootname
    #                                 WHERE headers.segment = 'FUVA' AND
    #                                     stims.stim1_x != -999 AND
    #                                     stims.stim1_y != -999 AND
    #                                     stims.stim2_x != -999 AND
    #                                     stims.stim2_y != -999;""")
    # data = [line for line in data]

    # y1 = [float(line.stim1_y) for line in data]
    # y2 = [float(line.stim2_y) for line in data]

    # im, nothin1, nothin2 = np.histogram2d(y2, y1, bins=200)
    # im = np.log(im)
    # ax.imshow(im, aspect='auto', interpolation='none')
    # ax.set_xlabel('y1')
    # ax.set_ylabel('y2')
    # ax.set_title('Segment A: Y vs Y')



    # ax = fig.add_subplot(2, 2, 3)
    # ax.grid(True)

    # data = engine.execute("""SELECT stim1_x, stim2_x
    #                                     FROM stims
    #                                     JOIN headers on stims.rootname = headers.rootname
    #                                     WHERE headers.segment = 'FUVB' AND
    #                                         stims.stim1_x != -999 AND
    #                                         stims.stim1_y != -999 AND
    #                                         stims.stim2_x != -999 AND
    #                                         stims.stim2_y != -999;""")
    # data = [line for line in data]

    # x1 = [float(line.stim1_x) for line in data]
    # x2 = [float(line.stim2_x) for line in data]

    # im, nothin1, nothin2 = np.histogram2d(x2, x1, bins=200)
    # im = np.log(im)
    # ax.imshow(im, aspect='auto', interpolation='none')
    # ax.set_xlabel('x1')
    # ax.set_ylabel('x2')
    # ax.set_title('Segment B: X vs X')



    # ax = fig.add_subplot(2, 2, 4)
    # ax.grid(True)

    # data = engine.execute("""SELECT stim1_y, stim2_y
    #                                     FROM stims
    #                                     JOIN headers on stims.rootname = headers.rootname
    #                                     WHERE headers.segment = 'FUVB' AND
    #                                         stims.stim1_x != -999 AND
    #                                         stims.stim1_y != -999 AND
    #                                         stims.stim2_x != -999 AND
    #                                         stims.stim2_y != -999;""")
    # data = [line for line in data]

    # y1 = [float(line.stim1_y) for line in data]
    # y2 = [float(line.stim2_y) for line in data]

    # im, nothin1, nothin2 = np.histogram2d(y2, y1, bins=200)
    # im = np.log(im)
    # ax.imshow(im, aspect='auto', interpolation='none')
    # ax.set_xlabel('y1')
    # ax.set_ylabel('y2')
    # ax.set_title('Segment B: Y vs Y')

    # #fig.colorbar(colors)
    # fig.savefig(os.path.join(out_dir, 'STIM_coord_relations_density.png'))
    # plt.close(fig)


    """
    print 1
    fig = plt.figure(1, figsize=(18, 12))
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.grid(True)
    x1 = [line[3] for line in data if '_a.fits' in line[0]]
    x2 = [line[5] for line in data if '_a.fits' in line[0]]
    times = [line[1] for line in data if '_a.fits' in line[0]]
    ax.scatter(x1, x2, times, s=5, c=times, alpha=.5, edgecolors='none')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_title('Segment A: X vs X')
    print 2
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    ax.grid(True)
    y1 = [line[4] for line in data if '_a.fits' in line[0]]
    y2 = [line[6] for line in data if '_a.fits' in line[0]]
    times = [line[1] for line in data if '_a.fits' in line[0]]
    ax.scatter(y1, y2, times, s=5, c=times, alpha=.5, edgecolors='none')
    slope, intercept = trend(y1, y2)
    print slope, intercept
    #plt.plot( [min(y1), max(y1)], [slope*min(y1)+intercept, slope*max(y1)+intercept], 'y--', lw=3)
    ax.set_xlabel('y1')
    ax.set_ylabel('y2')
    ax.set_title('Segment A: Y vs Y')
    print 3
    ax = fig.add_subplot(2, 2, 3, projection='3d')
    ax.grid(True)
    x1 = [line[3] for line in data if '_b.fits' in line[0]]
    x2 = [line[5] for line in data if '_b.fits' in line[0]]
    times = [line[1] for line in data if '_b.fits' in line[0]]
    ax.scatter(x1, x2, times, s=5, c=times, alpha=.5, edgecolors='none')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_title('Segment B: X vs X')
    print 4
    ax = fig.add_subplot(2, 2, 4, projection='3d')
    ax.grid(True)
    y1 = [line[4] for line in data if '_b.fits' in line[0]]
    y2 = [line[6] for line in data if '_b.fits' in line[0]]
    times = [line[1] for line in data if '_b.fits' in line[0]]
    colors = ax.scatter(y1, y2, times, s=5, c=times, alpha=.5, edgecolors='none')
    slope, intercept = trend(y1, y2)
    print slope, intercept
    #plt.plot( [min(y1), max(y1)], [slope*min(y1)+intercept, slope*max(y1)+intercept], 'y--', lw=3)
    ax.set_xlabel('y1')
    ax.set_ylabel('y2')
    ax.set_title('Segment B: Y vs Y')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar_ax.set_title('MJD')
    fig.colorbar(colors, cax=cbar_ax)
    fig.savefig(os.path.join(out_dir, 'STIM_coord_relations_time.png'))
    plt.close(fig)
    """

#-------------------------------------------------------------------------------
def send_email(missing_obs, missing_dates):
    """inform the parties that retrieval and calibration are done

    """
    sorted_index = np.argsort(np.array(missing_dates))
    missing_obs = np.array(missing_obs)[sorted_index]
    missing_dates = missing_dates[sorted_index]

    date_now = datetime.date.today()
    date_now = Time('%d-%d-%d' % (date_now.year, date_now.month, date_now.day), scale='utc', format='iso')
    date_diff = (date_now - missing_dates).sec / (60. * 60. * 24)

    index = np.where(date_diff < 7)[0]
    message = '--- WARNING ---\n'
    message += 'The following observations had missing stim pulses within the past week:\n'
    message += '-' * 40 + '\n'
    message += ''.join(['{} {}\n'.format(obs, date) for obs, date in zip(missing_obs[index], missing_dates.iso[index])])

    message += '\n\n\n'
    message += 'The following is a complete list of observations with missing stims.\n'
    message += '-' * 40 + '\n'
    message += ''.join(['{} {}\n'.format(obs, date) for obs, date in zip(missing_obs, missing_dates.iso)])

    svr_addr = 'smtp.stsci.edu'
    from_addr = 'ely@stsci.edu'
    #recipients = ['ely@stsci.edu', 'sahnow@stsci.edu', 'penton@stsci.edu', 'sonnentr@stsci.edu']
    recipients = 'mfix@stsci.edu'
    #to_addr = ', '.join(recipients)

    msg = MIMEMultipart()
    msg['Subject'] = 'Stim Monitor Report'
    msg['From'] = from_addr
    #msg['To'] = to_addr
    msg['To'] = recipients

    msg.attach(MIMEText(message))
    try:
        s = smtplib.SMTP(svr_addr)
        s.sendmail(from_addr, recipients, msg.as_string())
        s.quit()
    except:
        print("Cannot send email - perhaps you're not on the internal network?")

#-------------------------------------------------------------------------------

def stim_monitor():
    """Main function to monitor the stim pulses in COS observations

    1: populate the database
    2: find any datasets with missing stims [send email]
    3: make plots
    4: move plots to webpage
    5: check over individual observations

    """

    logger.info("STARTING MONITOR")

    settings = get_settings()

    webpage_dir = os.path.join(settings['webpage_location'], 'stim')
    monitor_dir = os.path.join(settings['monitor_location'], 'Stims')

    # for place in [webpage_dir, monitor_dir]:
    #     if not os.path.exists(place):
    #         logger.debug("creating monitor location: {}".format(place))
    #         os.makedirs(place)

    # missing_obs, missing_dates = find_missing()
    # send_email(missing_obs, missing_dates)

    # check_individual(monitor_dir, settings['connection_string'])
    # make_plots(monitor_dir, settings['connection_string'])
    
    logger.info("MAKING INTERACTIVE PLOTS")
    interactive_plotting(os.environ['HOME'], 'stim_peewee.html')
    # move_to_web(monitor_dir, webpage_dir)
    logger.info("FINISH MONITOR")