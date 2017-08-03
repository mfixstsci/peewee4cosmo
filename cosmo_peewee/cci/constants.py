"""Constants for monitor.py
"""

from __future__ import absolute_import
from datetime import datetime

__author__ = 'Justin Ely, Mees Fix'
__maintainer__ = 'Mees Fix'
__email__ = 'mfix@stsci.edu'
__status__ = 'Active'

__all__ = ['X_UNBINNED',
           'Y_UNBINNED',
           'X_BINNING',
           'Y_BINNING',
           'XLEN',
           'YLEN',
           'X_VALID',
           'Y_VALID',
           'FUVA_string',
           'FUVB_string',
           'MODAL_GAIN_LIMIT',
           'TIMESTAMP']

X_UNBINNED = 16384
Y_UNBINNED = 1024

#-- Optimal pore size for COS 'super pixel'..
X_BINNING = 8
Y_BINNING = 2  

XLEN = X_UNBINNED // X_BINNING
YLEN = Y_UNBINNED // Y_BINNING

X_VALID = (0 // X_BINNING, 16384 // X_BINNING)
Y_VALID = (0 // Y_BINNING, 1024 // Y_BINNING)

FUVA_string = '_00_'
FUVB_string = '_01_'

MODAL_GAIN_LIMIT = 3

date_time = str(datetime.now())
TIMESTAMP = (date_time.split()[0]+'T'+date_time.split()[1] ).replace(':','-')
