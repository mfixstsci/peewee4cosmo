from __future__ import absolute_import

import logging
logger = logging.getLogger(__name__)

from astropy.io import fits
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import numpy as np

import collections 
import itertools
import datetime
import seaborn as sns

from ..database.models import get_database, get_settings, Files
from ..database.models import Flagged_Pixels, Gain
from .constants import *

from bokeh.io import output_file, show, save
from bokeh.plotting import figure, ColumnDataSource
from bokeh.palettes import Category10
from bokeh.models import Span

from ..utils import remove_if_there

#-------------------------------------------------------------------------------
def make_overplot(gainsag_table, bm_hvlvl_a=167, bm_hvlvl_b=175, blue_modes=False):
    """
    Make plot of total gainmaps with gsagtab over plotted.

    Parameters
    ----------
    gainsag_table: str
        Path and filename of gainsag table
    bm_hvlvl_a: int
        HV level for FUVA bluemodes
    bm_hvlvl_b: int
        HV level for FUVB bluemodes
    blue_modes: bool
        Logic to make plots of blue modes or not.
    """

    #-- Get settings
    settings = get_settings()

    #-- Open GSAGTAB and gainmap.
    gsagtab = fits.open(gainsag_table)
    
    if blue_modes:
        fuva_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 'CCI', '*_00_{}*gainmap*'.format(bm_hvlvl_a)))
        fuva_gainmaps.sort()
        fuvb_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 'CCI', '*_01_{}*gainmap*'.format(bm_hvlvl_b)))
        fuvb_gainmaps.sort()
        
        filename = 'gainmap_gsag_overplot_bluemodes_{}.png'.format(datetime.date.today())

    else:
        #-- Sort all of the FUV gainmaps by segment
        fuva_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 'CCI', '*_00_*gainmap*'))
        fuva_gainmaps.sort()
        fuvb_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 'CCI', '*_01_*gainmap*'))
        fuvb_gainmaps.sort()
        
        filename = 'gainmap_gsag_overplot_{}.png'.format(datetime.date.today()) 

    #-- Open the last created gainmaps for each segment.
    fuva_gainmap = fits.open(fuva_gainmaps[-1:][0])
    fuvb_gainmap = fits.open(fuvb_gainmaps[-1:][0])

    #-- Set the HVLVL and SEGMENT from the A/B gainmaps
    hv_lvl_a = int(fuva_gainmap[0].header['DETHV'])
    hv_lvl_b = int(fuvb_gainmap[0].header['DETHV'])
    
    #-- Use the HV from the last A/B gainmaps to open the total gainmaps.
    fuva_total_gainmap = fits.open(os.path.join(settings['monitor_location'], 'CCI','total_gain_{}.fits'.format(hv_lvl_a)))
    fuvb_total_gainmap = fits.open(os.path.join(settings['monitor_location'], 'CCI','total_gain_{}.fits'.format(hv_lvl_b)))

    #-- gsagtab keywords for high voltage are HVLEVEL[A/B].
    hvlvl_key = {'FUVA':'HVLEVELA',
                 'FUVB':'HVLEVELB'}

    #-- Find which extension in the GSAGTAB matches the HVLVL and SEGMENT for the gainmap.
    for ext in range(1, len(gsagtab)):
        gsagtab_segment = gsagtab[ext].header['segment']
        if (hv_lvl_a == gsagtab[ext].header[hvlvl_key[gsagtab_segment]]) and (gsagtab[ext].header['segment']=='FUVA'):
            fuva_gsag_ext = ext
        elif (hv_lvl_b == gsagtab[ext].header[hvlvl_key[gsagtab_segment]]) and (gsagtab[ext].header['segment']=='FUVB'):
            fuvb_gsag_ext = ext
        else:
            pass

    #-- Create figure.
    f, axarr = plt.subplots(2, figsize=(20,15))

    #-- Set some plotting peramimeters.
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', lw=2)

    axes_font = 18
    title_font = 15

    #-- FUVA PANAL
    gm = axarr[0].imshow(fuva_total_gainmap['FUVALAST'].data, aspect='auto', cmap='gist_gray')
    f.colorbar(gm, ax=axarr[0])
    axarr[0].scatter(gsagtab[fuva_gsag_ext].data['LX'], gsagtab[fuva_gsag_ext].data['LY'], marker='+', color='red', label='DQ 8192')
    axarr[0].axhline(y=400, lw=2, ls='--', c='red')
    axarr[0].axhline(y=435, lw=2, ls='--', c='red')
    
    axarr[0].legend(fontsize=15)

    axarr[0].set_title('SEGMENT: FUVA, HVLEVEL: {}, Gainmap: {}'.format(hv_lvl_a, os.path.basename(fuva_gainmaps[-1:][0])), fontsize=title_font, fontweight='bold')
    axarr[0].set_xlabel('X (Pixels)', fontsize=axes_font, fontweight='bold')
    axarr[0].set_ylabel('Y (Pixels)', fontsize=axes_font, fontweight='bold')

    axarr[0].set_xlim([400, 15500])
    axarr[0].set_ylim([200 , 800])
    #-- END FUVA PANEL

    #-- FUVB PANEL
    gm = axarr[1].imshow(fuvb_total_gainmap['FUVBLAST'].data, aspect='auto', cmap='gist_gray')
    f.colorbar(gm, ax=axarr[1])
    axarr[1].scatter(gsagtab[fuvb_gsag_ext].data['LX'], gsagtab[fuvb_gsag_ext].data['LY'], marker='+', color='red', label='DQ 8192')
    axarr[1].axhline(y=460, lw=2, ls='--', c='red')
    axarr[1].axhline(y=495, lw=2, ls='--', c='red')
    
    axarr[1].legend(fontsize=15)
    
    axarr[1].set_title('SEGMENT: FUVB, HVLEVEL: {}, Gainmap: {}'.format(hv_lvl_b, os.path.basename(fuvb_gainmaps[-1:][0])), fontsize=title_font, fontweight='bold')
    axarr[1].set_xlabel('X (Pixels)', fontsize=axes_font, fontweight='bold')
    axarr[1].set_ylabel('Y (Pixels)', fontsize=axes_font, fontweight='bold')

    axarr[1].set_xlim([400, 15400])
    axarr[1].set_ylim([350, 800])
    #-- END FUVB PANEL

    #-- Save figure.
    plt.savefig(os.path.join(settings['monitor_location'], 'CCI', filename))
    #-- Close.
    plt.close()

#-------------------------------------------------------------------------------
def hotspot_plotter_interactive(segment):
    """Locate a plot hotspot's gain as a function of time.

    Parameters
    ----------
    segment: str
        FUVA or FUVB

    Returns
    -------
    None
    """

    #-- Get settings and connect to DB.
    settings = get_settings()
    database = get_database()
    database.connect()

    #-- HV histories for FUVA/FUVB
    hv_dictionary = {'FUVA':[163,167,169,171,173,178],
                     'FUVB':[163,167,169,175]}
    
    lp4_profile ={'FUVA':[400//Y_BINNING, 435//Y_BINNING],
                  'FUVB':[460//Y_BINNING, 495//Y_BINNING]}

    #-- Hot spot for FUVB B1
    hotspots = Gain.select(Gain.x, Gain.y).distinct().where(
                                                           (Gain.y.between(lp4_profile[segment][0], lp4_profile[segment][1])) &
                                                           (Gain.segment == segment) &
                                                           (Gain.gain <= 3)
                                                           )

    database.close()
    #-- Organize all x + y positions.
    coords = [(row.x, row.y) for row in hotspots]

    pix_history = list(Gain.select().where(
                                          (Gain.x.in_([row.x for row in hotspots])) &
                                          (Gain.y.in_([row.y for row in hotspots])) &
                                          (Gain.segment == segment)
                                          ).dicts())
    
    database.close()

    result = collections.defaultdict(list)
    
    #-- Organize all of the pixels by segment and .
    for d in pix_history:
        result[d['segment'], d['x'], d['y']].append(d)
    
    #-- If plots exist, delete to rewrite.
    plots = glob.glob(os.path.join(settings['monitor_location'],'CCI','hotspot_plots') + '/*html')
    if len(plots):
        for plot in plots:
            remove_if_there(plot)
            
    for combo in result.keys():
        segment,x,y = combo 
        filename = os.path.join(settings['monitor_location'],'CCI','hotspot_plots','hotspot_pixel_{}_{}_{}.html'.format(x*X_BINNING,y*Y_BINNING,segment))

        #-- Set up all the colors to cycle through for different HVs
        #-- https://bokeh.pydata.org/en/latest/docs/reference/palettes.html#bokeh-palettes
        palette = itertools.cycle(Category10[len(hv_dictionary[segment])])

        #-- Plotting config.
        output_file(filename)

        #-- Set plot height & width.
        plt_hgt = 1000
        plt_wth = 1500

        #-- Create figure
        p = figure(width=plt_wth, height=plt_hgt, title='({},{}) Gain VS Time.'.format(x*X_BINNING,y*Y_BINNING))
        p.title.text_font_size = '15pt'
        p.xaxis.axis_label_text_font_size = '15pt'
        p.yaxis.axis_label_text_font_size = '15pt'

        #-- Add lifetime positions and gainsag line.
        lp2 = Span(location=56101,
                    dimension='height', line_color='black',
                    line_dash='dashed', line_width=2)
        p.add_layout(lp2)

        lp3 = Span(location=57062,
                    dimension='height', line_color='black',
                    line_dash='dashed', line_width=2)
        p.add_layout(lp3)

        lp4 = Span(location=58028,
                    dimension='height', line_color='black',
                    line_dash='dashed', line_width=2)
        p.add_layout(lp4)

        sag_line = Span(location=3,
                    dimension='width', line_color='red',
                    line_dash='dashed', line_width=2)
        p.add_layout(sag_line)

        #-- For each HV, plot the gain vs time.
        for hv in hv_dictionary[segment]:
            gain = []
            date = []
            #-- Get gains and dates.   
            for dictionary in result[combo]:
                if dictionary['hv_lvl'] == hv:
                    gain.append(dictionary['gain'])
                    date.append(dictionary['expstart'])

            #-- Set titles and x+y labels
            p.xaxis.axis_label = 'Date (MJD)'
            p.yaxis.axis_label = "Modal Gain"
            p.circle(date, gain, legend=str(hv), size=6, color=next(palette), alpha=0.5)
        save(p, filename=filename)
#-------------------------------------------------------------------------------
