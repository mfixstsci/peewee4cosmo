"""Monitor HV Levels for COS FUV Detector
"""

import os

import matplotlib.pyplot as plt
import numpy as np

import scipy
from scipy.stats import linregress

from ..database.models import get_database, get_settings
from ..database.models import Hv_Level
from ..osm.monitor import fit_data

def plot_hv_levels(data):
    """Make HV vs time and HV residual checks for commanded HV and read back 
    HV.

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
    hv_lvl_in_volts = [d['hvlevel']* -15.69 - 2500 for d in data]
    

    # Plot HV vs Time
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15)
    plt.rc('axes', lw=2)

    f, (ax1, ax2, ax3) = plt.subplots(3, figsize=(20,10), sharex=True)
    f.suptitle('TIME VS. HVLVL | SEGMENT: {}'.format(segment[0]), fontsize=25)
    ax1.scatter(date, hv_lvl, c='k', facecolors='none', label='HV Step',
                alpha=0.5)
    ax1.set_ylabel('HV LEVEL (Steps)', fontsize=15, fontweight='bold')
    
    ax2.scatter(date, commanded_hv, c='k', facecolors='none', label='DETHVC',
                alpha=0.5)
    ax2.set_ylabel('HV LEVEL (Volts)', fontsize=15, fontweight='bold')
    
    ax3.scatter(date, read_back_hv, c='k', facecolors='none', label='DETHVL',
                alpha=0.5)
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

    # Plot HV Residuals vs time
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15)
    plt.rc('axes', lw=2)

    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(25,20))
    
    f.suptitle('HV RESIDUAL CHECK | SEGMENT: {}'.format(segment[0]), fontsize=25)
    ax1.scatter(date, hv_lvl_in_volts, c='k', facecolors='none', label='HV Step (Volts)',
                alpha=0.5)
    ax1.set_xlabel('Date (MJD)',fontsize=15,fontweight='bold')
    ax1.set_ylabel('HV STEP * -15.69 -2500 (Volts)', fontsize=15,fontweight='bold')
    ax2.scatter(date, commanded_hv, c='k', facecolors='none', label='DETHVC',
                alpha=0.5)
    ax2.set_xlabel('Date (MJD)',fontsize=15,fontweight='bold')
    ax2.set_ylabel('HV Commanded (Volts)', fontsize=15,fontweight='bold')
    fit,xdata,parameters,err = fit_data(hv_lvl_in_volts,
                                        commanded_hv)
    ax3.scatter(hv_lvl_in_volts, commanded_hv, c='k', facecolors='none',
                alpha=0.5)
    ax3.set_xlabel('HV STEP * -15.69 -2500 (Volts)',fontsize=15,fontweight='bold')
    ax3.set_ylabel('HV (Volts)', fontsize=15,fontweight='bold')
    
    ax3.plot(xdata, fit, c='r', ls='--', 
             label='Slope={}'.format(str(parameters[0])))

    ax4.scatter(date, np.array(commanded_hv) - np.array(hv_lvl_in_volts), c='k', facecolors='none', label='Residuals',
                alpha=0.5)
    ax4.set_xlabel('Date (MJD)',fontsize=15,fontweight='bold')
    ax4.set_ylabel('Residuals (Volts)', fontsize=15,fontweight='bold')
    for val in [-15.69, 15.69]:
        ax4.axhline(val, c='r', ls='-.')
    for ax in [ax1, ax2, ax3, ax4]:
        ax.grid(ls='-.')
        ax.legend(fontsize=15)

    figname = os.path.join(settings['monitor_location'],
                           'hvlvl',
                           'hv_in_volts_residuals_{}.png'.format(segment[0]))
    plt.savefig(figname)
    plt.close(f)
    


def main():
    database = get_database()
    database.connect()

    for segment in ['FUVA', 'FUVB']:

        data = Hv_Level.select().where(
                                       (Hv_Level.segment == segment)
                                     & (Hv_Level.hvlevel != -1)).dicts()
        plot_hv_levels(data)