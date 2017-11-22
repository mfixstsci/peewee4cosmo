"""
Set of utils that can access the DB through the command line. 
"""

import argparse
from peewee import *

from ..database.models import get_database, get_settings, Files
from ..database.models import Observations
from ..database.models import FUVA_raw_headers, FUVB_raw_headers, FUVA_corr_headers, FUVB_corr_headers
from ..database.models import NUV_corr_headers
#-------------------------------------------------------------------------------
def gather_all_data_rootname():
    """
    Provided a rootname, this returns all the data from all tables where rootname
    exists.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    settings = get_settings()
    database = get_database()
    
    database.connect()
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--rootname',
                        type=str,
                        help="Enter COS rootname")

    args = parser.parse_args()

    rootname = args.rootname

    detector = Observations.select(Observations.detector).where(Observations.rootname==rootname)

    for item in detector:
        detector = item.detector
    
    if detector == 'NUV':
        data = list(
                    (Observations.select()
                    .join(NUV_corr_headers, on=(Observations.rootname == NUV_corr_headers.rootname))
                    ).dicts())
    
    elif detector == 'FUV':
        data = list(
                    (Observations.select()
                    .join(FUVA_raw_headers, on=(Observations.rootname == FUVA_raw_headers.rootname))
                    .join(FUVB_raw_headers, on=(Observations.rootname == FUVB_raw_headers.rootname))
                    .join(FUVA_corr_headers, on=(Observations.rootname == FUVA_corr_headers.rootname))
                    .join(FUVB_corr_headers, on=(Observations.rootname == FUVB_corr_headers.rootname))
                    ).dicts())

    entry_to_print = filter(lambda dictionary: dictionary['rootname'] == rootname, data)
    database.close()
    print(entry_to_print)
#-------------------------------------------------------------------------------