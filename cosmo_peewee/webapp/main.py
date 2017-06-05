from flask import Flask, render_template, url_for

from peewee import *
import os
from ..database.models import get_settings, get_database
from ..database.models import Files
from ..database.models import Darks, Lampflash, Rawacqs

from ..dark.interactive_plots import plot_time

from ..osm.monitor import make_shift_table
from ..osm.monitor import make_interactive_plots as osm_interactive_plots
from ..dark.interactive_plots import plot_time as dark_interactive_plots
from ..dark.monitor import mjd_to_decyear 

from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8

from astropy.io import ascii

import numpy as np
#-------------------------------------------------------------------------------

app = Flask(__name__)

@app.route('/')
def index():
    return render_template("home.html")
#-------------------------------------------------------------------------------
@app.route('/search')
def search():
    return render_template("query.html")
#-------------------------------------------------------------------------------
@app.route('/monitoring')
def monitors():
   
    plots = []
    settings = get_settings()
    monitor_dir = settings['monitor_location']
    #-----------------------------------
    #-- OSM Interactive Plots
    plots = []
    settings = get_settings()

    osm_monitor_dir = os.path.join(settings['monitor_location'], 'Shifts')

    flash_data = make_shift_table(Lampflash)
    rawacq_data = make_shift_table(Rawacqs)

    #--grab the static resources
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    for detector in ['FUV', 'NUV']:
        fig = osm_interactive_plots(flash_data, rawacq_data, osm_monitor_dir, detector)
        #--render template
        script, div = components(fig)
        plots.append((script, div))
    #-----------------------------------
    #-- Dark Interactive Plots
    
    dark_monitor_dir = os.path.join(settings['monitor_location'], 'Darks')
    
    try:
        solar_data = np.genfromtxt(os.path.join(dark_monitor_dir, 'solar_flux.txt'), dtype=None)
        solar_date = np.array(mjd_to_decyear([line[0] for line in solar_data]))
        solar_flux = np.array([line[1] for line in solar_data])
    except TypeError:
        logger.warning("COULDN'T READ SOLAR DATA. PUTTING IN ZEROS.")
        solar_date = np.ones(1000)
        solar_flux = np.ones(1000)

    for detector in ['FUV', 'NUV']:
        if detector == 'FUV':
            segments = ['FUVA', 'FUVB']
        elif detector == 'NUV':
            segments = ['NUV']
        else:
            raise ValueError('Only FUV or NUV allowed. NOT:{}'.format(detector))
        
        for segment in segments:
            data = Darks.select().where(Darks.detector == segment)
                
            data = [row for row in data]
            dark = np.array([item.dark for item in data])
            mjd = np.array([item.date for item in data])
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
            outname = os.path.join(dark_monitor_dir,'dark_vs_time_{}.html'.format(segment))

            fig = dark_interactive_plots(detector, dark, mjd, temp, solar_flux, solar_date, outname)
            script, div = components(fig)
            plots.append((script, div))
    #-----------------------------------
    #-- Make plots appear in website.
    html = render_template(
        'monitors.html',
        plots=plots,
        js_resources=js_resources,
        css_resources=css_resources,
    )
    return encode_utf8(html)
#-------------------------------------------------------------------------------
def run():
    app.run()
#-------------------------------------------------------------------------------
def run_debug():
    app.run(debug=True)
