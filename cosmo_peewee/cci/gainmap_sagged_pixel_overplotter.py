from __future__ import absolute_import

"""Gainmap/GSAGTAB over plotting + plotting script.
"""

import os

import argparse
from astropy.io import fits
from astropy.time import Time
from bokeh.io import output_file, show, save
from bokeh.models import Span
from bokeh.plotting import figure, ColumnDataSource
from bokeh.palettes import Category10
import collections 
import datetime
import functools
import glob
import itertools
import logging
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import multiprocessing as mp
import numpy as np
import seaborn as sns

from .constants import *
from ..database.models import get_database, get_settings, Files
from ..database.models import Flagged_Pixels, Gain
from ..utils import remove_if_there

logger = logging.getLogger(__name__)


def make_overplot(gainsag_table, bm_hvlvl_a=167, bm_hvlvl_b=175, blue_modes=False):
    """Make plot of total gainmaps with gsagtab over plotted.

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

    # Get settings
    settings = get_settings()

    # Open GSAGTAB and gainmap.
    gsagtab = fits.open(gainsag_table)
    
    if blue_modes:
        fuva_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 
                                               'CCI', '*_00_{}*gainmap*'\
                                                    .format(bm_hvlvl_a)))
        fuva_gainmaps.sort()
        fuvb_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 
                                               'CCI', '*_01_{}*gainmap*'\
                                                    .format(bm_hvlvl_b)))
        fuvb_gainmaps.sort()
        
        filename = 'gainmap_gsag_overplot_bluemodes_{}.png'\
                        .format(datetime.date.today())

    else:
        # Sort all of the FUV gainmaps by segment
        fuva_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 
                                               'CCI', '*_00_*gainmap*'))
        fuva_gainmaps.sort()
        fuvb_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 
                                               'CCI', '*_01_*gainmap*'))
        fuvb_gainmaps.sort()
        
        filename = 'gainmap_gsag_overplot_{}.png'.format(datetime.date.today())

    # Open the last created gainmaps for each segment.
    fuva_gainmap = fits.open(fuva_gainmaps[-1:][0])
    fuvb_gainmap = fits.open(fuvb_gainmaps[-1:][0])

    # Set the HVLVL and SEGMENT from the A/B gainmaps
    hv_lvl_a = int(fuva_gainmap[0].header['DETHV'])
    hv_lvl_b = int(fuvb_gainmap[0].header['DETHV'])
    
    # Use the HV from the last A/B gainmaps to open the total gainmaps.
    fuva_total_gainmap = fits.open(os.path.join(settings['monitor_location'], 
                                                'CCI','total_gain_{}.fits'\
                                                    .format(hv_lvl_a)))
    fuvb_total_gainmap = fits.open(os.path.join(settings['monitor_location'], 
                                                'CCI','total_gain_{}.fits'\
                                                    .format(hv_lvl_b)))

    # gsagtab keywords for high voltage are HVLEVEL[A/B].
    hvlvl_key = {'FUVA':'HVLEVELA',
                 'FUVB':'HVLEVELB'}

    # Find which extension in the GSAGTAB matches the HVLVL and SEGMENT for 
    # the gainmap.
    for ext in range(1, len(gsagtab)):
        gsagtab_segment = gsagtab[ext].header['segment']
        if (hv_lvl_a == gsagtab[ext].header[hvlvl_key[gsagtab_segment]]) and \
                        (gsagtab[ext].header['segment']=='FUVA'):
            fuva_gsag_ext = ext
        elif (hv_lvl_b == gsagtab[ext].header[hvlvl_key[gsagtab_segment]]) and \
                        (gsagtab[ext].header['segment']=='FUVB'):
            fuvb_gsag_ext = ext
        else:
            pass

    # Create figure.
    f, axarr = plt.subplots(2, figsize=(20,15))

    # Set some plotting peramimeters.
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', lw=2)

    axes_font = 18
    title_font = 15

    # FUVA PANAL
    gm = axarr[0].imshow(fuva_total_gainmap['FUVALAST'].data, aspect='auto', 
                         cmap='gist_gray')
    f.colorbar(gm, ax=axarr[0])
    axarr[0].scatter(gsagtab[fuva_gsag_ext].data['LX'], 
                     gsagtab[fuva_gsag_ext].data['LY'], 
                     marker='+', color='red', 
                     label='DQ 8192')
    axarr[0].axhline(y=400, lw=2, ls='--', c='red')
    axarr[0].axhline(y=435, lw=2, ls='--', c='red')
    
    axarr[0].legend(fontsize=15)

    axarr[0].set_title('SEGMENT: FUVA, HVLEVEL: {}, Gainmap: {}'\
                            .format(hv_lvl_a, os.path.basename(
                                                fuva_gainmaps[-1:][0])), 
                       fontsize=title_font, fontweight='bold')
    axarr[0].set_xlabel('X (Pixels)', fontsize=axes_font, fontweight='bold')
    axarr[0].set_ylabel('Y (Pixels)', fontsize=axes_font, fontweight='bold')

    axarr[0].set_xlim([400, 15500])
    axarr[0].set_ylim([380 ,480])
    # END FUVA PANEL

    # FUVB PANEL
    gm = axarr[1].imshow(fuvb_total_gainmap['FUVBLAST'].data, aspect='auto', cmap='gist_gray')
    f.colorbar(gm, ax=axarr[1])
    axarr[1].scatter(gsagtab[fuvb_gsag_ext].data['LX'], 
                     gsagtab[fuvb_gsag_ext].data['LY'], 
                     marker='+', color='red', 
                     label='DQ 8192')
    axarr[1].axhline(y=460, lw=2, ls='--', c='red')
    axarr[1].axhline(y=495, lw=2, ls='--', c='red')
    
    axarr[1].legend(fontsize=15)
    
    axarr[1].set_title('SEGMENT: FUVB, HVLEVEL: {}, Gainmap: {}'\
                       .format(hv_lvl_b, 
                       os.path.basename(fuvb_gainmaps[-1:][0])), 
                       fontsize=title_font, fontweight='bold')
    axarr[1].set_xlabel('X (Pixels)', fontsize=axes_font, fontweight='bold')
    axarr[1].set_ylabel('Y (Pixels)', fontsize=axes_font, fontweight='bold')

    axarr[1].set_xlim([400, 15400])
    axarr[1].set_ylim([450, 530])
    # END FUVB PANEL

    # Save figure.
    plt.savefig(os.path.join(settings['monitor_location'], 'CCI', filename))
    # Close.
    plt.close()


def gsagtab_overplot_comparison(hv_lvl, compare=False, potential_gsagtab=None, current_gsagtab=None):
    """Compare two gain sag tables to see gain sag progression.
    Convention is to have latest gain sag table be plotted with red + symbols 
    and the older table be plotted in green + symbols. This can accept any two 
    gsagtabs but user must make sure that the potential_gsagtab is created after
    current_gsagtab to ensure proper output.

    Parameters
    ----------
    potential_gsagtab: str
        path to gsagtab you want to use in pipeline
    current_gsagtab: str
        path to gsagtab currently in use
    """

    if compare:
        filename = 'gsagtab_comparison_{}.pdf'.format(hv_lvl)
    else:
        filename = 'gainmap_gsagtab_{}_overplot.pdf'.format(hv_lvl)

    # Give blue modes different filename to not overwrite.    
    if '_blue.fits' in os.path.basename(potential_gsagtab):
        filename = filename.replace('.pdf', '_blue.pdf')

    settings = get_settings()
    
    # Open gsagtabs
    potential_hdu = fits.open(potential_gsagtab)
    if compare:
        current_hdu = fits.open(current_gsagtab)

    # Check the hdu lengths to make sure they are the same size.
    if compare:
        assert (len(potential_hdu) == len(current_hdu)), \
                'GSAGTABS ARE DIFFERENT LENGTH'
    
    gainmap = fits.open(os.path.join(settings['monitor_location'], 
                                     'CCI','total_gain_{}.fits'.format(hv_lvl)))

    # gsagtab keywords for high voltage are HVLEVEL[A/B].
    hvlvl_key = {'FUVA':'HVLEVELA',
                 'FUVB':'HVLEVELB'}

    # Find which extension in the GSAGTAB matches the HVLVL and SEGMENT for the gainmap.
    for ext in range(1, len(potential_hdu)):
        gsagtab_segment = potential_hdu[ext].header['segment']
        if (hv_lvl == potential_hdu[ext].header[hvlvl_key[gsagtab_segment]]) and \
           (potential_hdu[ext].header['segment']=='FUVA'):
            fuva_gsag_ext = ext
        elif (hv_lvl == potential_hdu[ext].header[hvlvl_key[gsagtab_segment]]) and \
             (potential_hdu[ext].header['segment']=='FUVB'):
            fuvb_gsag_ext = ext
        else:
            pass
    

    # Create figure.
    f, axarr = plt.subplots(2, figsize=(20,15))

    # Set some plotting peramimeters.
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', lw=2)

    axes_font = 18
    title_font = 15


    if compare:
        legend_title = 'Latest'
    else:
        legend_title = 'DQ 8192'

    # An x1d for LP4 with the 2-zone extraction performed.
    # This is a 1222 exposure, which has the widest profile. 
    # I will pull the encircled enegry contours from this figure.
    profiles = fits.open(os.path.join(settings['monitor_location'],'CCI', 
                                            'ldel01p6q_x1d.fits'))

    # Find respective segment information.
    profile_FUVA = np.where(profiles[1].data['segment'] == 'FUVA')[0][0]
    profile_FUVB = np.where(profiles[1].data['segment'] == 'FUVB')[0][0]

    # FUVA PANAL
    gm = axarr[0].imshow(gainmap['FUVALAST'].data, aspect='auto', cmap='gist_gray')
    f.colorbar(gm, ax=axarr[0])
    
    axarr[0].scatter(potential_hdu[fuva_gsag_ext].data['LX'], 
                     potential_hdu[fuva_gsag_ext].data['LY'], 
                     marker='+', color='red', 
                     label=legend_title)
    if compare:
        axarr[0].scatter(current_hdu[fuva_gsag_ext].data['LX'], 
                         current_hdu[fuva_gsag_ext].data['LY'], 
                         marker='+', color='green', 
                         label='Current')

    # List for header keywords to plot contours
    contours = ['y_upper_outer', 'y_lower_outer', 'y_lower_inner', 
                'y_upper_inner']
    contour_colors = ['r','r','g','g']
    
    #-- Loop through and plot the contours with their respective colors.
    for contour,color in zip(contours, contour_colors):
        axarr[0].plot(np.arange(len(profiles[1].data[contour][profile_FUVA])), 
                      profiles[1].data[contour][profile_FUVA], c=color)
    
    # axarr[0].axhline(y=400, lw=2, ls='--', c='red')
    # axarr[0].axhline(y=435, lw=2, ls='--', c='red')
    
    axarr[0].legend(fontsize=15)

    axarr[0].set_title('SEGMENT: FUVA, HVLEVEL: {}, Gainmap: {}'\
                            .format(hv_lvl, 'total_gain_{}.fits'.format(hv_lvl)), 
                        fontsize=title_font, fontweight='bold')
    axarr[0].set_xlabel('X (Pixels)', fontsize=axes_font, fontweight='bold')
    axarr[0].set_ylabel('Y (Pixels)', fontsize=axes_font, fontweight='bold')

    axarr[0].set_xlim([400, 15500])
    axarr[0].set_ylim([200 , 800])
    # END FUVA PANEL

    # FUVB PANEL
    gm = axarr[1].imshow(gainmap['FUVBLAST'].data, aspect='auto', 
                         cmap='gist_gray')
    f.colorbar(gm, ax=axarr[1])
    
    axarr[1].scatter(potential_hdu[fuvb_gsag_ext].data['LX'],
                     potential_hdu[fuvb_gsag_ext].data['LY'], 
                     marker='+', color='red', 
                     label=legend_title)
    if compare:
        axarr[1].scatter(current_hdu[fuvb_gsag_ext].data['LX'], 
                         current_hdu[fuvb_gsag_ext].data['LY'], 
                         marker='+', color='green', 
                         label='Current')
    
    for contour,color in zip(contours, contour_colors):
        axarr[1].plot(np.arange(len(profiles[1].data[contour][profile_FUVB])), 
                      profiles[1].data[contour][profile_FUVB], c=color)

    # axarr[1].axhline(y=460, lw=2, ls='--', c='red')
    # axarr[1].axhline(y=495, lw=2, ls='--', c='red')
    
    axarr[1].legend(fontsize=15)
    
    axarr[1].set_title('SEGMENT: FUVB, HVLEVEL: {}, Gainmap: {}'\
                            .format(hv_lvl, 'total_gain_{}.fits'.format(hv_lvl)), 
                       fontsize=title_font, fontweight='bold')
    axarr[1].set_xlabel('X (Pixels)', fontsize=axes_font, fontweight='bold')
    axarr[1].set_ylabel('Y (Pixels)', fontsize=axes_font, fontweight='bold')

    axarr[1].set_xlim([400, 15400])
    axarr[1].set_ylim([350, 800])
    # END FUVB PANEL

    # Plot enclosed energy contours for the blue modes at LP2.
    # FUVA is at hv_lvl 173
    if hv_lvl == 173:
        # 1096
        # blue_profiles = fits.open(os.path.join(settings['monitor_location'],'CCI', 'lddy01hqq_x1d.fits.gz'))
        # 1055 
        blue_profiles = fits.open(os.path.join(settings['monitor_location'],
                                               'CCI', 
                                               'ldcv20esq_x1d.fits.gz'))
        # Find respective segment information.
        blue_profile_FUVA = np.where(blue_profiles[1].data['segment'] \
                                        == 'FUVA')[0][0]        
        for contour,color in zip(contours, contour_colors):
            axarr[0].plot(np.arange(len(blue_profiles[1].data[contour][blue_profile_FUVA])), 
                          blue_profiles[1].data[contour][blue_profile_FUVA], c=color)

    # FUVB is at hv_lvl 175
    if hv_lvl == 175:
        # 1096
        # blue_profiles = fits.open(os.path.join(settings['monitor_location'],'CCI', 'lddy01hqq_x1d.fits.gz'))
        # 1055 
        blue_profiles = fits.open(os.path.join(settings['monitor_location'],
                                               'CCI', 
                                               'ldcv20esq_x1d.fits.gz'))
        # Find respective segment information.
        blue_profile_FUVB = np.where(blue_profiles[1].data['segment'] \
                                        == 'FUVB')[0][0]
        for contour,color in zip(contours, contour_colors):
            axarr[1].plot(np.arange(len(blue_profiles[1].data[contour][blue_profile_FUVB])),
                          blue_profiles[1].data[contour][blue_profile_FUVB], c=color)


    # Save figure, change permissions.
    comparison_filename = os.path.join(settings['monitor_location'], 'CCI',
                             'gsagtab_comparisons', filename)
    remove_if_there(comparison_filename)
    plt.savefig(comparison_filename)
    os.chmod(comparison_filename, 0o776)

    # Close.
    plt.close()


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

    settings = get_settings()
    database = get_database()
    database.connect()

    # HV histories for FUVA/FUVB
    hv_dictionary = {'FUVA':[163,167,169,171,173,178],
                     'FUVB':[163,167,169,175]}
    
    lp4_profile ={'FUVA':[400//Y_BINNING, 440//Y_BINNING],
                  'FUVB':[460//Y_BINNING, 500//Y_BINNING]}

    # Hot spot for FUVB B1
    hotspots = Gain.select(Gain.x, Gain.y)\
                    .distinct().where(
                        (Gain.y.between(lp4_profile[segment][0], 
                                        lp4_profile[segment][1]))
                        & (Gain.segment == segment)
                        & (Gain.gain <= 3))

    database.close()
    # Organize all x + y positions.
    coords = [(row.x, row.y) for row in hotspots]

    pix_history = list(Gain.select()
                            .where(
                                (Gain.x.in_([row.x for row in hotspots]))
                                & (Gain.y.in_([row.y for row in hotspots]))
                                & (Gain.segment == segment)).dicts())
    
    database.close()

    result = collections.defaultdict(list)
    
    # Organize all of the pixels by segment and .
    for d in pix_history:
        result[d['segment'], d['x'], d['y']].append(d)
    
    # If plots exist, delete to rewrite.
    plots = glob.glob(os.path.join(settings['monitor_location'],
                                   'CCI','hotspot_plots') 
                                   + '/*{}.html'.format(segment))
    if len(plots):
        for plot in plots:
            remove_if_there(plot)
            
    for combo in result.keys():
        segment,x,y = combo 
        filename = os.path.join(settings['monitor_location'],
                                'CCI',
                                'hotspot_plots','hotspot_pixel_{}_{}_{}.html'\
                                    .format(x*X_BINNING,y*Y_BINNING,segment))

        # Set up all the colors to cycle through for different HVs
        # https://bokeh.pydata.org/en/latest/docs/reference/palettes.html#bokeh-palettes
        palette = itertools.cycle(Category10[len(hv_dictionary[segment])])

        # Plotting config.
        output_file(filename)
        os.chmod(filename, 0o776)

        # Set plot height & width.
        plt_hgt = 1000
        plt_wth = 1500

        # Create figure
        p = figure(width=plt_wth, height=plt_hgt, 
                   title='({},{}) Gain VS Time.'\
                    .format(x*X_BINNING,y*Y_BINNING))
        p.title.text_font_size = '15pt'
        p.xaxis.axis_label_text_font_size = '15pt'
        p.yaxis.axis_label_text_font_size = '15pt'

        # Add lifetime positions and gainsag line.
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

        # For each HV, plot the gain vs time.
        for hv in hv_dictionary[segment]:
            gain = []
            date = []
            # Get gains and dates.   
            for dictionary in result[combo]:
                if dictionary['hv_lvl'] == hv:
                    gain.append(dictionary['gain'])
                    date.append(dictionary['expstart'])

            # Set titles and x+y labels
            p.xaxis.axis_label = 'Date (MJD)'
            p.yaxis.axis_label = "Modal Gain"
            p.circle(date, gain, legend=str(hv), size=6, 
                     color=next(palette), alpha=0.5)
        save(p, filename=filename)


def gsagtab_plot_by_date():
    """Plot gain sag holes onto gainmaps.

    Parameters:
    -----------
    segment: str
        FUVA or FUVB
    min_date: float
        Lower range for dates to search for
    max_date:
        Upper range for dates to search for
    
    Returns
    -------
    None
    """

    settings = get_settings()
    database = get_database()
    
    # Current Time
    now = Time(datetime.datetime.now().isoformat()).mjd
        
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--segment',
                        type=str,
                        default='FUVB',
                        help="FUVA or FUVB")
    
    #-- For dates, defaults give sagging activity over the past 10 days.
    parser.add_argument('--min_date',
                        type=float,
                        default=now-14.0,
                        help="Minimum date when gainsag holes appeared")
    
    parser.add_argument('--max_date',
                        type=float,
                        default=now,
                        help="Minimum date when gainsag holes appeared")
    
    parser.add_argument('--compare',
                        type=bool,
                        default=False,
                        help="Compare to CRDS gsagtab")

    args = parser.parse_args()
    
    sagged_pix = list(Flagged_Pixels.select()
                            .where(
                                (Flagged_Pixels.mjd_below_3.between(args.min_date,
                                                                    args.max_date))
                                & (Flagged_Pixels.segment == args.segment)).dicts())
    database.close()
    result = collections.defaultdict(list)
    
    # Organize all of the pixels by hv_lvls.
    for d in sagged_pix:
        result[d['hv_lvl']].append(d)
    
    hvlvl_key = {'FUVA':'HVLEVELA',
                 'FUVB':'HVLEVELB'}
    
    for hv in result.keys():
        if hv < 163:
            continue
        
        if args.compare:
            crds_gsagtab = fits.open(os.path.join(settings['lref'], 
                                     'zbn1927gl_gsag.fits'))

            for ext in range(1,len(crds_gsagtab)):
                try:
                    if (crds_gsagtab[ext].header[hvlvl_key[args.segment]] == hv) and \
                       (crds_gsagtab[ext].header['segment']==args.segment):
                        lx = crds_gsagtab[ext].data['LX']
                        ly = crds_gsagtab[ext].data['LY']
                    else:
                        continue
                except KeyError as e:
                    print(e)
        
        if args.compare:
            lx = [dictionary['x']*X_BINNING for dictionary in sagged_pix]
            ly = [dictionary['y']*Y_BINNING for dictionary in sagged_pix]
        
        else:
            x = [dictionary['x']*X_BINNING for dictionary in result[hv]]
            y = [dictionary['y']*Y_BINNING for dictionary in result[hv]]
        
        gainmap = fits.open(os.path.join(settings['monitor_location'],'CCI','total_gain_{}.fits'.format(hv)))

        plt.figure(figsize=(20,10))
        plt.rc('xtick', labelsize=15) 
        plt.rc('ytick', labelsize=15)
        plt.rc('axes', lw=2)    
        
        plt.imshow(gainmap['{}LAST'.format(args.segment)].data, aspect='auto', cmap='gist_gray')
        plt.scatter(x, y, marker='+', color='red', label='HV {}'.format(hv))
        
        if args.compare:
            plt.scatter(lx, ly, marker='+', color='green', label='CRDS HV {}'.format(hv))
        
        plt.title('Segment {} Gain Map {} For Dates {} -- {}'.format(args.segment,hv, args.min_date, args.max_date), fontsize=15)
        
        if args.segment == 'FUVA':
            plt.xlim(400, 15500)
            plt.ylim(200 , 800)
        if args.segment == 'FUVB':
            plt.xlim(400, 15400)
            plt.ylim(350, 800)
        
        plt.xlabel('XCORR (Pixels)', fontsize=20)
        plt.ylabel('YCORR (Pixels)', fontsize=20)
        plt.legend(fontsize=15)
        
        if args.compare:
            filename='gsag_by_date_compare_{}-{}_{}_{}.png'.format(int(args.min_date), int(args.max_date), hv, args.segment)
        else:
            filename='gsag_by_date_{}-{}_{}_{}.png'.format(int(args.min_date), int(args.max_date), hv, args.segment)
        
        plt.savefig(os.path.join(settings['monitor_location'], 'CCI', 'gsagtab_comparisons', filename))
        plt.close()
        print(os.path.join(settings['monitor_location'], 'CCI', 'gsagtab_comparisons', filename))


def compare_and_plot_gsagtable_data(ext, gsagtab_old=None, gsagtab_new=None, outdir=None):
    """Compare data from two gain sag tables and plot differences.

    Parameters
    ----------
    gsagtab_old: str
        Path to gainsag table with eariler date
    gsagtab_new: str
        Path to gainsag table that has date > gsagtab_old
    """

    # settings = get_settings()

    segment_key = {'FUVA':'FUVALAST',
                   'FUVB':'FUVBLAST'}
    
    hv_lvl_key = {'FUVA':'HVLEVELA',
                  'FUVB':'HVLEVELB'}

    hdu_old = fits.open(gsagtab_old)
    hdu_new = fits.open(gsagtab_new)

    
    segment = hdu_new[ext].header['segment']
    hv_lvl = hdu_new[ext].header[hv_lvl_key[segment]]
    
    if hv_lvl not in [163,167,169,171,173,175,178]:
        return
    
    if segment == 'FUVA':
        xlim = [400, 15500]
        ylim = [200 , 800]
    else:
        xlim = [400, 15400]
        ylim = [350, 800]
    
    # Set some plotting peramimeters.
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', lw=2)

    f, (ax1, ax2) = plt.subplots(2, figsize=(25,15))
    
    # Find coordinates that aren't in 
    new_coords = [(x,y) for x,y in zip(hdu_new[ext].data['LX'], hdu_new[ext].data['LY'])]
    old_coords = [(x,y) for x,y in zip(hdu_old[ext].data['LX'], hdu_old[ext].data['LY'])]
    
    gainmap = fits.open(os.path.join('/grp/hst/cos/Monitors/', 'CCI','total_gain_{}.fits'.format(hv_lvl)))
    
    ax1.imshow(gainmap[segment_key[segment]].data, aspect='auto', cmap='gist_gray')
    ax1.scatter(hdu_new[ext].data['LX'], hdu_new[ext].data['LY'], marker='+', c='r', label=os.path.basename(gsagtab_new))
    ax1.scatter(hdu_old[ext].data['LX'], hdu_old[ext].data['LY'], marker='x', c='g', label=os.path.basename(gsagtab_old))
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_xlabel('X (Pixels)', fontsize=18, fontweight='bold')
    ax1.set_ylabel('Y (Pixels)', fontsize=18, fontweight='bold')
    ax1.set_title('SEGMENT: {} | HVLEVEL: {}'.format(segment, hv_lvl), fontsize=20, fontweight='bold')
    ax1.legend(fontsize=13)
    
    ax2.imshow(gainmap[segment_key[segment]].data, aspect='auto', cmap='gist_gray')
    
    residuals = list(set(new_coords) - set(old_coords))
    if residuals:
        diff_x, diff_y = zip(*residuals)
        print(segment, hv_lvl, zip([x/8 for x in diff_x], [y/2 for y in diff_y]))
        ax2.scatter(diff_x, diff_y, marker='+', c='m', label='Difference')
        ax2.legend(fontsize=13)

    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    ax2.set_xlabel('X (Pixels)', fontsize=18, fontweight='bold')
    ax2.set_ylabel('Y (Pixels)', fontsize=18, fontweight='bold')
    ax2.set_title('SEGMENT: {} | HVLEVEL: {}'.format(segment, hv_lvl), fontsize=20, fontweight='bold')
            

    filename = 'gsagtab_residual_comparion_{}_{}.png'.format(segment, hv_lvl)
    plt.savefig(os.path.join(outdir, filename))
    plt.close()


def compare_and_plot_gsagtable_data_entry():
    """
    Entry Point for plots comparing two gainsag tables.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    settings = get_settings()
    pool = mp.Pool(processes=settings['num_cpu'])

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--old_gsagtab',
                        type=str,
                        help="Path to gsagtab")
    
    parser.add_argument('--new_gsagtab',
                        type=str,
                        help="Path to gsagtab")
    
    parser.add_argument('--out_dir',
                        type=str,
                        help="Path you want to write plot out to.")
    
    args = parser.parse_args()

    partial = functools.partial(compare_and_plot_gsagtable_data,
                                gsagtab_old=args.old_gsagtab, 
                                gsagtab_new=args.new_gsagtab, 
                                outdir=args.out_dir)

    pool.map(partial, range(1,78))


def plot_gainmap_and_gsagtab_by_hv(gsagtab, hv_lvl):
    """Give a HV and gain sag tables, plot segments FUVA + FUVB in single plot.

    Parameters
    ----------
    gsagtab: str
        Path to gsagtab
    hv_lvl: int
        HV level you are interested in.
    """

    settings = get_settings()
    out_dir = os.path.join(settings['monitor_location'], 'CCI', 'gainmap_gsagtab_delivery_plots')
    
    gsagtab_hdu = fits.open(gsagtab)
    gainmap = fits.open(os.path.join(settings['monitor_location'], 'CCI', 'total_gain_{}.fits'.format(hv_lvl)))
    
    # gsagtab keywords for high voltage are HVLEVEL[A/B].
    hvlvl_key = {'FUVA':'HVLEVELA',
                 'FUVB':'HVLEVELB'}

    # Find which extension in the GSAGTAB matches the HVLVL and SEGMENT for the gainmap.
    for ext in range(1, len(gsagtab_hdu)):
        gsagtab_segment = gsagtab_hdu[ext].header['segment']
        if (hv_lvl == gsagtab_hdu[ext].header[hvlvl_key[gsagtab_segment]]) and \
            (gsagtab_hdu[ext].header['segment']=='FUVA'):
            fuva_gsag_ext = ext
        elif (hv_lvl == gsagtab_hdu[ext].header[hvlvl_key[gsagtab_segment]]) and \
            (gsagtab_hdu[ext].header['segment']=='FUVB'):
            fuvb_gsag_ext = ext
        else:
            pass
    
    # Set some plotting peramimeters.
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', lw=2)
    
    gsagtab_name = os.path.basename(gsagtab)

    filename = 'gainmap_{}_{}.png'.format(os.path.splitext(gsagtab_name)[0],hv_lvl)
    f, (ax1, ax2) = plt.subplots(2, figsize=(25,15))
    
    plt.subplots_adjust(hspace=0.5)

    
    f.suptitle('GSAGTAB: {}'.format(gsagtab_name), fontsize=25, fontweight='bold')
    
    gm = ax1.imshow(gainmap['FUVALAST'].data, aspect='auto', cmap='gist_gray')
    ax1.scatter(gsagtab_hdu[fuva_gsag_ext].data['LX'], 
                gsagtab_hdu[fuva_gsag_ext].data['LY'], 
                marker='+', c='r', label='FUVA {}'.format(hv_lvl))
    ax1.set_title('FUVA GAINMAP + GSAG', fontsize=20, fontweight='bold')
    ax1.set_xlabel('XCORR (Pixels)', fontsize=20, fontweight='bold')
    ax1.set_ylabel('YCORR (Pixels)', fontsize=20, fontweight='bold')
    ax1.set_xlim([0,16384])
    ax1.set_ylim([300,600])
    ax1.legend(fontsize=20)
    f.colorbar(gm, ax=ax1)

    gm = ax2.imshow(gainmap['FUVBLAST'].data, aspect='auto', cmap='gist_gray')
    ax2.scatter(gsagtab_hdu[fuvb_gsag_ext].data['LX'], 
                gsagtab_hdu[fuvb_gsag_ext].data['LY'], 
                marker='+', c='r', label='FUVB {}'.format(hv_lvl))
    ax2.set_title('FUVB GAINMAP + GSAG', fontsize=20, fontweight='bold')
    ax2.set_xlabel('XCORR (Pixels)', fontsize=20, fontweight='bold')
    ax2.set_ylabel('YCORR (Pixels)', fontsize=20, fontweight='bold')
    ax2.set_xlim([0,16384])
    ax2.set_ylim([300,800])
    ax2.legend(fontsize=20)
    f.colorbar(gm, ax=ax2)
    
    plt.savefig(os.path.join(out_dir, filename))
    plt.close(f)


def plot_gainmap_and_gsagtab_by_hv_entry():
    """
    Entry Point for plot_gainmap_and_gsagtab_by_hv.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    parser = argparse.ArgumentParser()    
    parser.add_argument('--gsagtab',
                        type=str,
                        help="Path to gsagtab")
    
    parser.add_argument('--hv_lvl',
                        type=int,
                        help="High Voltage Level")
    args = parser.parse_args()

    plot_gainmap_and_gsagtab_by_hv(args.gsagtab, args.hv_lvl)
