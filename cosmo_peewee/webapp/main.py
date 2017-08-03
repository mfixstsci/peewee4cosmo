from peewee import *
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

import os

import numpy as np

from jinja2 import Markup
from flask import Flask, Response, redirect, json, render_template, request, session, url_for
from wtforms import Form, BooleanField, StringField, validators
from werkzeug.utils import secure_filename
import wtforms

import datetime

from form_options import FORM_OPTIONS
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
def select_multi_checkbox(field, ul_class='checky-boiz', **kwargs):
    """
    Function to allow for multiple checked boxes.

    Parameters
    ----------
    field: str
        Name of the checkbox.
    ul_class: str
        something....
    
    Returns
    -------
    wtforms object.
    """
    kwargs.setdefault('type', 'checkbox')
    field_id = kwargs.pop('id', field.id)
    html = [u'<ul %s>' % wtforms.widgets.html_params(id=field_id, class_=ul_class)]
    for value, label, checked in field.iter_choices():
        choice_id = u'%s-%s' % (field_id, value)
        options = dict(kwargs, name=field.name, value=value, id=choice_id)
        if checked:
            options['checked'] = 'checked'
        html.append(u'<li><input %s />&nbsp;' % wtforms.widgets.html_params(**options))
        html.append(u'<label for="%s">%s</label></li>' % (choice_id, label))
    html.append(u'</ul>')
    return Markup(u''.join(html))
#-------------------------------------------------------------------------------
app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret key'
@app.route('/')
def index():
    """
    Returns the 'home.html template'

    Parameters
    ----------
    None

    Returns
    -------
    home.html webpage.
    """
    return render_template("home.html")
#-------------------------------------------------------------------------------

class QueryForm(wtforms.Form):
    """
    Creates WTForm....

    Parameters
    ----------
    wtforms.Form: WTForms Form Object

    Returns
    -------
    None
    """
    rootname = wtforms.StringField()
    #-- Detectors
    detectors = wtforms.SelectMultipleField(choices=FORM_OPTIONS['detector'], 
                                            widget=select_multi_checkbox)
    #-- Observation Types
    obstypes = wtforms.SelectMultipleField(choices=FORM_OPTIONS['obstype'], 
                                           widget=select_multi_checkbox)
    #-- Imgage Types
    img_type = wtforms.SelectMultipleField(choices=FORM_OPTIONS['imagetyp'], 
                                           widget=select_multi_checkbox)
    #-- Gratings
    gratings = wtforms.SelectMultipleField(choices=FORM_OPTIONS['grating'], 
                                           widget=select_multi_checkbox)
    #-- Cenwaves
    cenwaves = wtforms.SelectMultipleField(choices=FORM_OPTIONS['cenwave'], 
                                           widget=select_multi_checkbox)
    #-- Lifetime Adjustments
    life_adj = wtforms.SelectMultipleField(choices=FORM_OPTIONS['life_adj'], 
                                           widget=select_multi_checkbox)
    #-- Exposure Types
    exptype = wtforms.SelectMultipleField(choices=FORM_OPTIONS['exptype'], 
                                           widget=select_multi_checkbox)
#-------------------------------------------------------------------------------
@app.route('/search')
def search():
    """
    Build and return query page.

    Parameters
    ----------
    None

    Returns
    -------
    query page template.
    """
    form = QueryForm(request.args)
    if request.query_string and form.validate(): # validate if any arguments were submitted

        #-- Parse all fields and only select selected fields from form.
        #-- form.data returns a dictionary
        print('Flask is fun.....')
    return render_template("query.html", form=form)
#-------------------------------------------------------------------------------
@app.route('/monitoring')
def monitors():
    """
    Returns the monitoring page

    Parameters
    ----------
    None

    Returns
    -------
    Monitor template page.
    """
    return render_template('monitors.html')
    
#-------------------------------------------------------------------------------
@app.route('/monitoring/lampflash_monitor')
def osm_monitor():
    """
    Renders the Lampflash monitoring plots

    Parameters
    ----------
    None

    Returns
    -------
    Lampflash monitors template.
    """
    
    plots = []
    settings = get_settings()
    #monitor_dir = settings['monitor_location']
    
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
    """
    Renders the dark monitor plots.

    Parameters
    ----------
    None

    Returns
    -------
    Dark monitor plots.
    """

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

            outname = os.path.join(dark_monitor_dir, 'dark_vs_time_{}.html'.format(segment))

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
    """
    Run Flask server

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    app.run()
#-------------------------------------------------------------------------------
def run_debug():
    """
    Run Flask server in debug mode for testing and development. 
    
    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    app.run(debug=True)
