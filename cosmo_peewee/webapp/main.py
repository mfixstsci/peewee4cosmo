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

from flask import Flask, Response, redirect, json, render_template, request, session, url_for
from wtforms import Form, BooleanField, StringField, validators
from werkzeug.utils import secure_filename
import wtforms
from flask_wtf import FlaskForm
from flask_wtf.file import FileField as FlaskFileField, FileRequired as FlaskFileRequired

import datetime

#-------------------------------------------------------------------------------
#-- Global Constants
DETECTORS = (
    ('fuv', 'FUV'),
    ('nuv', 'NUV'),
)

OBS_TYPE = (
    ('imaging', 'Imaging'),
    ('spectroscopic', 'Spectroscopic'),
)

IMG_TYPE = (
    ('timetag', 'Time-Tag'),
    ('accum', 'Accum'),
)
#-------------------------------------------------------------------------------
def select_multi_checkbox(field, ul_class='', **kwargs):
    kwargs.setdefault('type', 'checkbox')
    field_id = kwargs.pop('id', field.id)
    html = [u'<ul class="checky-boiz" %s>' % wtforms.widgets.html_params(id=field_id, class_=ul_class)]
    for value, label, checked in field.iter_choices():
        choice_id = u'%s-%s' % (field_id, value)
        options = dict(kwargs, name=field.name, value=value, id=choice_id)
        if checked:
            options['checked'] = 'checked'
        html.append(u'<li><input %s /> ' % wtforms.widgets.html_params(**options))
        html.append(u'<label for="%s">%s</label></li>' % (choice_id, label))
    html.append(u'</ul>')
    return u''.join(html)
#-------------------------------------------------------------------------------
app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret key'
@app.route('/')
def index():
    return render_template("home.html")
#-------------------------------------------------------------------------------

class ExampleForm(FlaskForm):
    rootname = wtforms.StringField()
    #-- Detectors
    # fuv = wtforms.BooleanField()
    # nuv = wtforms.BooleanField()
    detectors = wtforms.SelectMultipleField(choices=DETECTORS, widget=select_multi_checkbox)
    #-- Obs Type
    obs_type = wtforms.SelectMultipleField(choices=OBS_TYPE, widget=select_multi_checkbox)
    #-- Img Type
    img_type = wtforms.SelectMultipleField(choices=IMG_TYPE, widget=select_multi_checkbox)
    #-- Gratings
    g130m = wtforms.BooleanField() # FUV
    g140l = wtforms.BooleanField()
    g160m = wtforms.BooleanField()
    g185m = wtforms.BooleanField() # NUV
    g225m = wtforms.BooleanField()
    g285m = wtforms.BooleanField()
    g230l = wtforms.BooleanField()
    #-- Proposal ID field.
    proposid = wtforms.StringField()
    apertures = wtforms.widgets.Select()
    #-- Cenwaves FUV
    cw_1055 = wtforms.BooleanField()# G130M
    cw_1096 = wtforms.BooleanField()
    cw_1222 = wtforms.BooleanField()
    cw_1291 = wtforms.BooleanField()
    cw_1300 = wtforms.BooleanField()
    cw_1309 = wtforms.BooleanField()
    cw_1318 = wtforms.BooleanField()
    cw_1327 = wtforms.BooleanField()
    cw_1577 = wtforms.BooleanField()# G160M
    cw_1589 = wtforms.BooleanField()
    cw_1600 = wtforms.BooleanField()
    cw_1611 = wtforms.BooleanField()
    cw_1623 = wtforms.BooleanField()
    cw_1105 = wtforms.BooleanField()# G140L

#-------------------------------------------------------------------------------
@app.route('/search', methods=['GET', 'POST'])
def search():
    form = ExampleForm()  # FlaskForm knows how to reach into `request` for data, unlike wtforms.Form
    if request.method == 'GET':
        print('GET')
        
    return render_template("query.html",form=form)
#-------------------------------------------------------------------------------
@app.route('/monitoring')
def monitors():
    return render_template('monitors.html')
    
#-------------------------------------------------------------------------------
@app.route('/monitoring/lampflash_monitor')
def osm_monitor():
    plots = []
    settings = get_settings()
    monitor_dir = settings['monitor_location']
    
    #--grab the static resources
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    #-----------------------------------
    #-- OSM Interactive Plots
    osm_monitor_dir = os.path.join(settings['monitor_location'], 'Shifts')

    flash_data = make_shift_table(Lampflash)
    rawacq_data = make_shift_table(Rawacqs)

    for detector in ['FUV', 'NUV']:
        fig = osm_interactive_plots(flash_data, rawacq_data, osm_monitor_dir, detector)
        #--render template
        script, div = components(fig)
        plots.append((script, div))
    
    html = render_template(
        'individual_monitor.html',
        plots=plots,
        js_resources=js_resources,
        css_resources=css_resources,
        monitor_name='Lampflash Monitor',
    )
    return encode_utf8(html)
#-------------------------------------------------------------------------------
@app.route('/monitoring/dark_monitor')
def dark_monitor():
    
    plots = []
    settings = get_settings()
    monitor_dir = settings['monitor_location']

    #--grab the static resources
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()
    
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

            index_keep = np.where((longitude < 250) | (latitude > 10))[0]
            
            mjd = mjd[index_keep]
            dark = dark[index_keep]
            temp = temp[index_keep]
            latitude = latitude[index]
            longitude = longitude[index]

            outname = os.path.join(dark_monitor_dir,'dark_vs_time_{}.html'.format(segment))

            fig = dark_interactive_plots(detector, dark, mjd, temp, solar_flux, solar_date, outname)
            script, div = components(fig)
            plots.append((script, div))
    #-----------------------------------
    #-- Make plots appear in website.
    html = render_template(
        'individual_monitor.html',
        plots=plots,
        js_resources=js_resources,
        css_resources=css_resources,
        monitor_name='Dark Monitor',
    )
    return encode_utf8(html)
#-------------------------------------------------------------------------------
def run():
    app.run()
#-------------------------------------------------------------------------------
def run_debug():
    app.run(debug=True)
