from __future__ import absolute_import

"""Routine to monitor the modal gain in each pixel as a
function of time.  Uses COS Cumulative Image (CCI) files
to produce a modal gain map for each time period.  Modal gain
maps for each period are collated to monitor the progress of
each pixel(superpixel) with time.  Pixels that drop below
a threshold value are flagged and collected into a
gain sag table reference file (gsagtab).

The PHA modal gain threshold is set by global variable MODAL_GAIN_LIMIT.
Allowing the modal gain of a distribution to come within 1 gain bin
of the threshold results in ~8% loss of flux.  Within
2 gain bins, ~4%
3 gain bins, ~2%
4 gain bins, ~1%

However, due to the column summing, a 4% loss in a region does not appear to be so in the extracted spectrum.
"""

__author__ = 'Justin Ely, Mees Fix'
__maintainer__ = 'Mees Fix'
__email__ = 'mfix@stsci.edu'
__status__ = 'Active'

import argparse
import os
import shutil
import sys
import time
from datetime import datetime
import gzip
import glob
import logging
logger = logging.getLogger(__name__)
import time

from astropy.io import fits
from astropy.modeling import models, fitting
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import leastsq, newton, curve_fit
import fitsio

from ..utils import rebin, enlarge
from .constants import * 

if sys.version_info.major == 2:
    from itertools import izip as zip

from ..database.models import get_database, get_settings, Files

from copy import deepcopy
#------------------------------------------------------------

class CCI:
    """
    Creates a cci_object designed for use in the monitor.

    Each COS cumulative image fits file is made into its
    own cci_object where the data is more easily used.
    Takes a while to make, contains a few large arrays and
    various header keywords from the original file.
    """

    def __init__(self, filename, **kwargs):
        """Open filename and create CCI Object"""

        #-- Assign full filename
        self.input_file = filename

        #-- Set x and y bins for COS super pixels
        self.xbinning = kwargs.get('xbinning', 1)
        self.ybinning = kwargs.get('ybinning', 1)

        #-- For a super pixel, it must have atleast 30 counts to be considered for fitting.
        self.mincounts = kwargs.get('mincounts', 30)

        #-- Strip the path off of filename
        path, cci_name = os.path.split(filename)

        #-- Split the ext name off of the filename
        cci_name, ext = os.path.splitext(cci_name)

        #-- trim off any remaining extensions: .fits, .gz, .tar, etc
        while not ext == '':
            cci_name, ext = os.path.splitext(cci_name)

        self.cci_name = cci_name

        #-- Open the file.
        try:
        	self.open_fits()
        except:
            print('COULD NOT COMPLETE open_fits METHOD IN CCI FOR {}'.format(self.input_file))
            database = get_database()
            database.connect()

            query = Files.update(monitor_flag=0).where(Files.filename==os.path.basename(self.input_file))
            query.execute()

            database.close()

            return
         	
        #-- CCI's are comprised of ~ 1 week of COS exposure time
        #-- We want to make sure that the CCI actually has COS exposures in it.
        if not self.numfiles:
            print('numfiles IS EMPTY, NO DATA IN CCI {}'.format(self.input_file))
            database = get_database()
            database.connect()
            
            query = Files.update(monitor_flag=0).where(Files.filename==os.path.basename(self.input_file))
            query.execute()
            
            database.close()
            
            return

        #-- CCI's are 3 dimensional (x, y, PHA). For each x,y pixel position there are
        #-- 31 PHA extensions. Each the counts in each PHA extension are binned into a 
        #-- histogram that its then fit for each pixel. The mean of this fit is the 
        #-- modal gain value for that pixel. The gain at each pixel construct the gainmap.
        gainmap, counts, std = measure_gainimage(self.big_array)
        self.gain_image = gainmap
        self.std_image = std

        #-------------------------------
        # NOT SURE IF WE USE THIS BLOCK!
        #-------------------------------

        # if kwargs.get('only_active_area', True):
        #     brftab = os.path.join(os.environ['lref'], self.brftab)
        #     left, right, top, bottom = read_brftab(brftab, self.segment)

        #     left //= self.xbinning
        #     right //= self.xbinning

        #     top //= self.ybinning
        #     bottom //= self.ybinning

        #     self.gain_image[:bottom] = 0
        #     self.gain_image[top:] = 0
        #     self.gain_image[:, :left] = 0
        #     self.gain_image[:, right:] = 0

        # if kwargs.get('ignore_spots', True):
        #     ### Dynamic when delivered to CRDS
        #     reffiles = glob.glob(os.path.join(os.environ['lref'], '*spot.fits'))
        #     creation_dates = np.array([fits.getval(item, 'DATE') for item in reffiles])
        #     spottab = reffiles[creation_dates.argmax()]

        #     if os.path.exists(spottab):
        #         regions = read_spottab(spottab,
        #                                self.segment,
        #                                self.expstart,
        #                                self.expend)

        #         for lx, ly, dx, dy in regions:
        #             lx //= self.xbinning
        #             dx //= self.xbinning

        #             ly //= self.ybinning
        #             dy //= self.ybinning

        #             #-- pad the regions by 1 bin in either direction
        #             lx -= 1
        #             dx += 2
        #             ly -= 1
        #             dy += 2
        #             #--

        #             self.gain_image[ly:ly+dy, lx:lx+dx] = 0

        # self.gain_index = np.where(self.gain_image > 0)
        # self.bad_index = np.where((self.gain_image <= 3) &
        #                           (self.gain_image > 0))

    def open_fits(self):
        """Open CCI file and populated attributes with
        header keywords and data arrays.
        """

        #-- Open fits file.
        hdu = fits.open(self.input_file)
        primary = hdu[0].header

        #-- Make sure that the shape of CCI are correct.
        try:
            assert (hdu[2].data.shape == (Y_UNBINNED, X_UNBINNED)), 'ERROR: Input CCI not standard dimensions'
        except AssertionError, e:
            print(self.input_file)
            raise Exception(e.args)

        #-- Assign a bunch of attributes
        self.detector = primary['DETECTOR']
        self.segment = primary['SEGMENT']
        self.obsmode = primary['OBSMODE']
        self.expstart = primary['EXPSTART']
        self.expend = primary['EXPEND']
        self.exptime = primary['EXPTIME']
        self.numfiles = primary['NUMFILES']
        self.counts = primary['COUNTS']
        self.dethv = primary.get('DETHV', -999)
        try:
            self.brftab = primary['brftab'].split('$')[1]
        except:
            self.brftab = 'x1u1459il_brf.fits'

        if self.expstart:
            #-- Finds to most recently created HVTAB
            settings = get_settings()
            hvtable_list = glob.glob(os.path.join(settings['lref'], '*hv.fits'))
            HVTAB = hvtable_list[np.array([fits.getval(item, 'DATE') for item in hvtable_list]).argmax()]
            hvtab = fits.open(HVTAB)

            if self.segment == 'FUVA':
                hv_string = 'HVLEVELA'
            elif self.segment == 'FUVB':
                hv_string = 'HVLEVELB'

        self.file_list = [line[0].decode("utf-8") for line in hdu[1].data]
        self.big_array = np.array([rebin(hdu[i+2].data, bins=(self.ybinning, self.xbinning)) for i in range(32)])
        self.get_counts(self.big_array)
        self.extracted_charge = self.pha_to_coulombs(self.big_array)
        

        self.gain_image = np.zeros((YLEN, XLEN))
        self.modal_gain_width = np.zeros((YLEN, XLEN))

        self.cnt00_00 = len(self.big_array[0].nonzero()[0])
        self.cnt01_01 = len(self.big_array[1].nonzero()[0])
        self.cnt02_30 = len(self.big_array[2:31].nonzero()[0])
        self.cnt31_31 = len(self.big_array[31:].nonzero()[0])
        

    def get_counts(self, in_array):
        """collapse pha arrays to get total counts accross all
        PHA bins.

        Will also search for and add in accum data if any exists.
        """

        out_array = np.sum(in_array, axis=0)

        #-- Test before implementation
        #-- Should only effect counts and charge extensions.
        #-- no implications for modal gain arrays or measurements
        if self.segment == 'FUVA':
            accum_name = self.cci_name.replace('00_','02_')  ##change when using OPUS data
        elif self.segment == 'FUVB':
            accum_name = self.cci_name.replace('01_','03_')  ##change when using OPUS data
        else:
            accum_name = None
            print('ERROR: name not standard')

        if os.path.exists(accum_name):
            accum_data = rebin(fits.getdata(CCI_DIR+accum_name, 0),bins=(Y_BINNING,self.xbinning))
            out_array += accum_data
            self.accum_data = accum_data
        else:
            self.accum_data = None

        self.counts_image = out_array

    def pha_to_coulombs(self, in_array):
        """Convert pha to picocoloumbs to calculate extracted charge.
        Equation comes from D. Sahnow.
        """

        coulomb_value = 1.0e-12*10**((np.array(range(0,32))-11.75)/20.5)
        zlen, ylen, xlen = in_array.shape
        out_array = np.zeros((ylen, xlen))

        for pha,layer in enumerate(in_array):
            out_array += (coulomb_value[pha]*layer)

        return out_array

    def write(self, out_name=None):
        '''Write current CCI object to fits file.

        Output files are used in later analysis to determine when
        regions fall below the threshold.
        '''

        out_name = out_name or self.cci_name + '_gainmap.fits'

        if os.path.exists(out_name):
            print("{} EXISTS, NOT CLOBBERING".format(out_name))
            return

        #-------Ext=0
        hdu_out = fits.HDUList(fits.PrimaryHDU())

        hdu_out[0].header['TELESCOP'] = 'HST'
        hdu_out[0].header['INSTRUME'] = 'COS'
        hdu_out[0].header['DETECTOR'] = 'FUV'
        hdu_out[0].header['OPT_ELEM'] = 'ANY'
        hdu_out[0].header['FILETYPE'] = 'GAINMAP'

        hdu_out[0].header['XBINNING'] = self.xbinning
        hdu_out[0].header['YBINNING'] = self.ybinning
        hdu_out[0].header['SRC_FILE'] = self.cci_name
        hdu_out[0].header['SEGMENT'] = self.segment
        hdu_out[0].header['EXPSTART'] = self.expstart
        hdu_out[0].header['EXPEND'] = self.expend
        hdu_out[0].header['EXPTIME'] = self.exptime
        hdu_out[0].header['NUMFILES'] = self.numfiles
        hdu_out[0].header['COUNTS'] = self.counts
        hdu_out[0].header['DETHV'] = self.dethv
        hdu_out[0].header['cnt00_00'] = self.cnt00_00
        hdu_out[0].header['cnt01_01'] = self.cnt01_01
        hdu_out[0].header['cnt02_30'] = self.cnt02_30
        hdu_out[0].header['cnt31_31'] = self.cnt31_31

        #-------EXT=1
        included_files = np.array(self.file_list)
        files_col = fits.Column('files', '24A', 'rootname', array=included_files)
        tab = fits.BinTableHDU.from_columns([files_col])

        hdu_out.append(tab)
        hdu_out[1].header['EXTNAME'] = 'FILES'

        #-------EXT=2
        hdu_out.append(fits.ImageHDU(data=self.gain_image))
        hdu_out[2].header['EXTNAME'] = 'MOD_GAIN'

        #-------EXT=3
        hdu_out.append(fits.ImageHDU(data=self.counts_image))
        hdu_out[3].header['EXTNAME'] = 'COUNTS'

        #-------EXT=4
        hdu_out.append(fits.ImageHDU(data=self.extracted_charge))
        hdu_out[4].header['EXTNAME'] = 'CHARGE'

        #-------EXT=5
        hdu_out.append(fits.ImageHDU(data=self.big_array[0]))
        hdu_out[5].header['EXTNAME'] = 'cnt00_00'

        #-------EXT=6
        hdu_out.append(fits.ImageHDU(data=self.big_array[1]))
        hdu_out[6].header['EXTNAME'] = 'cnt01_01'

        #-------EXT=7
        hdu_out.append(fits.ImageHDU(data=np.sum(self.big_array[2:31],axis=0)))
        hdu_out[7].header['EXTNAME'] = 'cnt02_30'

        #-------EXT=8
        hdu_out.append(fits.ImageHDU(data=self.big_array[31]))
        hdu_out[8].header['EXTNAME'] = 'cnt31_31'

        #-------Write to file
        hdu_out.writeto(out_name)
        hdu_out.close()

#-------------------------------------------------------------------------------

def make_total_gainmap(hv_lvl, gainmap_dir=None, segment='FUV', start_mjd=55055, end_mjd=70000, reverse=False):

    #-- Depending on the segment, the filename are different.
    if segment == 'FUVA':
        search_string = 'l_*_00_{}_cci_gainmap.fits*'.format(hv_lvl)
    elif segment == 'FUVB':
        search_string = 'l_*_01_{}_cci_gainmap.fits*'.format(hv_lvl)

    #-- Get all of the data and sort.
    all_datasets = [item for item in glob.glob(os.path.join(gainmap_dir, search_string))]
    all_datasets.sort()
    
    #-- Switch the order to reverse to see what the detector originally looked like.
    if reverse:
        all_datasets = all_datasets[::-1]

    #-- Make an empty array the size of the detector.
    out_data = np.zeros((YLEN, XLEN))

    #-- For each gainmap
    for item in all_datasets:
        
        cci_hdu = fits.open(item)

        #-- Make sure you are making a total gainmap in the dates you defined.
        if not cci_hdu[0].header['EXPSTART'] >= start_mjd:
            continue
        if not cci_hdu[0].header['EXPSTART'] <= end_mjd:
            continue
        
        #-- Get the data
        cci_data = cci_hdu['MOD_GAIN'].data
        
        #-- Find where there is data 
        index = np.where(cci_data)
        
        #-- get the high voltage level.
        dethv = cci_hdu[0].header['DETHV']

        #-- Add the data to the array and move to the next gainmap.
        out_data[index] = cci_data[index]

    return enlarge(out_data, x=X_BINNING, y=Y_BINNING)

#-------------------------------------------------------------------------------
def make_all_gainmaps(gainmap_dir, start_mjd=55055, end_mjd=70000, total=False):
    """
    Make all of the total gainmaps.

    Parameters
    ----------
    gainmap_dir: str
        Directory of where the gainmaps live
    start_mjd: int
        Start date of for gainmaps you want to build. (Default will catch first gainmap)
    end_mjd: int
        End date that gainmaps should have. (Default will be date way after recent gainmap)
    total: bool
        Make gainmap from all HV levels. 
    """

    if total:
        filename = os.path.join(gainmap_dir,'total_gain.fits')
        hdu_out = fits.HDUList(fits.PrimaryHDU())

        #-- Adding primary header with file specifications to make results reproducible

        hdu_out[0].header['TELESCOP'] = 'HST'
        hdu_out[0].header['INSTRUME'] = 'COS'
        hdu_out[0].header['DETECTOR'] = 'FUV'
        hdu_out[0].header['OPT_ELEM'] = 'ANY'
        hdu_out[0].header['FILETYPE'] = 'GAINMAP'

        hdu_out[0].header['EXPSTART'] = start_mjd
        hdu_out[0].header['EXP_END'] = end_mjd
        hdu_out[0].header['CCI_DIR'] = gainmap_dir
        hdu_out[0].header['HVLEVEL'] = 'ALL'
        hdu_out[0].header['DATEMADE'] = (time.strftime("%d/%m/%Y"))

        #-- Data ext
        hdu_out.append(fits.ImageHDU(data=make_total_gainmap('???', gainmap_dir, 'FUVA', start_mjd, end_mjd, reverse=True)))
        hdu_out[1].header['EXTNAME'] = 'FUVAINIT'
        hdu_out.append(fits.ImageHDU(data=make_total_gainmap('???', gainmap_dir, 'FUVB', start_mjd, end_mjd, reverse=True)))
        hdu_out[2].header['EXTNAME'] = 'FUVBINIT'
        hdu_out.append(fits.ImageHDU(data=make_total_gainmap('???', gainmap_dir, 'FUVA', start_mjd, end_mjd)))
        hdu_out[3].header['EXTNAME'] = 'FUVALAST'
        hdu_out.append(fits.ImageHDU(data=make_total_gainmap('???', gainmap_dir, 'FUVB', start_mjd, end_mjd)))
        hdu_out[4].header['EXTNAME'] = 'FUVBLAST'
        hdu_out.writeto(filename, clobber=True)
        hdu_out.close()
    else:
        for hv in range(150, 179):
            filename = os.path.join(gainmap_dir,'total_gain_{}.fits'.format(hv))
            hdu_out = fits.HDUList(fits.PrimaryHDU())

            #-- Adding primary header with file specifications to make results reproducible

            hdu_out[0].header['TELESCOP'] = 'HST'
            hdu_out[0].header['INSTRUME'] = 'COS'
            hdu_out[0].header['DETECTOR'] = 'FUV'
            hdu_out[0].header['OPT_ELEM'] = 'ANY'
            hdu_out[0].header['FILETYPE'] = 'GAINMAP'

            hdu_out[0].header['EXPSTART'] = start_mjd
            hdu_out[0].header['EXP_END'] = end_mjd
            hdu_out[0].header['CCI_DIR'] = gainmap_dir
            hdu_out[0].header['HVLEVEL'] = hv
            hdu_out[0].header['DATEMADE'] = (time.strftime("%d/%m/%Y"))

            #-- Data ext
            hdu_out.append(fits.ImageHDU(data=make_total_gainmap(hv, gainmap_dir, 'FUVA', start_mjd, end_mjd, reverse=True)))
            hdu_out[1].header['EXTNAME'] = 'FUVAINIT'
            hdu_out.append(fits.ImageHDU(data=make_total_gainmap(hv, gainmap_dir, 'FUVB', start_mjd, end_mjd, reverse=True)))
            hdu_out[2].header['EXTNAME'] = 'FUVBINIT'
            hdu_out.append(fits.ImageHDU(data=make_total_gainmap(hv, gainmap_dir, 'FUVA', start_mjd, end_mjd)))
            hdu_out[3].header['EXTNAME'] = 'FUVALAST'
            hdu_out.append(fits.ImageHDU(data=make_total_gainmap(hv, gainmap_dir, 'FUVB', start_mjd, end_mjd)))
            hdu_out[4].header['EXTNAME'] = 'FUVBLAST'
            hdu_out.writeto(filename, clobber=True)
            hdu_out.close()

#-------------------------------------------------------------------------------

def write_and_pull_gainmap(cci_name, out_dir='None'):
    """Make modal gainmap for cos cumulative image.

    """

    #-- Get settings
    settings = get_settings()

    #-- Set output directory for gainmaps.
    out_dir = os.path.join(settings['monitor_location'],'CCI')
    
    #-- Set full file_path from peewee object and run CCI class to obtain gainmap
    full_path = os.path.join(cci_name.path, cci_name.filename)
    current = CCI(full_path, xbinning=X_BINNING, ybinning=Y_BINNING)
    out_name = os.path.join(out_dir, os.path.basename(full_path.replace('.fits.gz', '_gainmap.fits.gz')))

    logger.debug("WRITING GAINMAP TO: {}".format(out_name))
    current.write(out_name)

    index = np.where(current.gain_image > 0)

    info = {'filename': os.path.basename(full_path),
            'segment': current.segment,
            'hv_lvl': int(current.dethv),
            'expstart': round(current.expstart, 5)}
    
    #-- If the y indices are empty, x will be empty as well.
    #-- Just make an entry that is zeros.
    for y, x in zip(*index):
        info['gain'] = round(float(current.gain_image[y, x]), 3)
        info['counts'] = round(float(current.counts_image[y, x]), 3)
        info['std_dev'] = round(float(current.std_image[y, x]), 3)
        info['x'] = int(x)
        info['y'] = int(y)

        yield deepcopy(info)

#-------------------------------------------------------------------------------