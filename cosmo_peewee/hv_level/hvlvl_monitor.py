"""Monitor HV Levels for COS FUV Detector
"""

import os

import matplotlib.pyplot as plt
import numpy as np

import scipy
from scipy.stats import linregress

from ..database.models import get_database, get_settings
from ..database.models import Hv_Level


def plot_hv_levels(data):
    """Make different HV plots

    Parameters
    ----------
    data: list
        List of dictionaries containing data

    Returns
    -------
    None
    """
    
    settings = get_settings()

    date = [d['expstart'] for d in data]
    segment = [d['segment'] for d in data]
    read_back_hv = [d['dethvl'] for d in data]
    commanded_hv = [d['dethvc'] for d in data]
    hv_lvl = [d['hvlevel'] for d in data]

    # Plot commanded vs read back.
    f, (ax1, ax2) = plt.subplots(2, figsize=(20,10), sharex=True)
    f.suptitle('COMMAND HV VS READ BACK HV | SEGMENT: {}'.format(segment[0]), fontsize=25)
    ax1.scatter(read_back_hv, commanded_hv, c='k')
    fit,xdata,parameters,err = fit_data(read_back_hv,
                                        commanded_hv)
    ax1.plot(xdata, fit, c='r', ls='--', 
             label='Slope={}'.format(str(parameters[0])))
    ax1.set_ylabel('COMMANDED HV (VOLTS)', fontsize=15, fontweight='bold')
    ax1.legend(fontsize=15)
    residuals = commanded_hv - fit
    ax2.scatter(xdata, residuals, c='k', label='Residuals')
    ax2.set_xlabel('READ BACK HV (VOLTS)', fontsize=15, fontweight='bold')
    ax2.set_ylabel('COMMANDED HV - Fit (Volts)', fontsize=15, fontweight='bold')
    for val in [-15.89, 15.89]:
        ax2.axhline(val, c='r', ls='-.')
    ax2.legend(fontsize=15)
    figname = os.path.join(settings['monitor_location'],
                           'hvlvl',
                           'command_vs_readback_{}.png'.format(segment[0]))
    plt.savefig(figname)
    plt.close(f)

    # Plot HV vs Time
    f, (ax1, ax2, ax3) = plt.subplots(3, figsize=(20,10), sharex=True)
    f.suptitle('TIME VS. HVLVL | SEGMENT: {}'.format(segment[0]), fontsize=25)
    ax1.scatter(date, hv_lvl, c='k', facecolors='none', label='HV Step')
    ax1.set_ylabel('HV LEVEL (Steps)', fontsize=15, fontweight='bold')
    ax2.scatter(date, commanded_hv, c='k', facecolors='none', label='DETHVC')
    ax2.set_ylabel('HV LEVEL (Volts)', fontsize=15, fontweight='bold')
    ax3.scatter(date, read_back_hv, c='k', facecolors='none', label='DETHVL')
    ax3.set_xlabel('Date (MJD)', fontsize=20, fontweight='bold')
    ax3.set_ylabel('HV LEVEL (Volts)', fontsize=15, fontweight='bold')
    for ax in [ax1, ax2, ax3]:
        ax.grid(ls='-.')
        ax.legend(fontsize=15)
    figname = os.path.join(settings['monitor_location'], 
                           'hvlvl', 
                           'hv_vs_time_{}.png'.format(segment[0]))
    plt.savefig(figname)
    plt.close(f)


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


def main():
    database = get_database()
    database.connect()

    for segment in ['FUVA', 'FUVB']:

        data = Hv_Level.select().where(
                                       (Hv_Level.segment == segment)
                                     & (Hv_Level.hvlevel != -1)).dicts()
        plot_hv_levels(data)