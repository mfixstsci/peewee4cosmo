"""Set of utils that can access the DB through the command line. 
"""

import argparse
import os
import pandas as pd
from peewee import *


from ..database.models import get_database, get_settings, Files
from ..database.models import Observations
from ..database.models import FUVA_raw_headers, FUVB_raw_headers, FUVA_corr_headers, FUVB_corr_headers
from ..database.models import NUV_corr_headers



def gather_all_data_rootname():
    """
    Provided a rootname, this returns all the data from all tables where 
    rootname exists.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    database = get_database()
    
    database.connect()
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--rootname',
                        type=str,
                        help="Enter COS rootname")

    args = parser.parse_args()

    rootname = args.rootname

    detector = Observations.select(Observations.detector).where(
                                   Observations.rootname == rootname)

    for item in detector:
        detector = item.detector
    
    if detector == 'NUV':
        data = list(
                    (Observations.select()
                    .join(NUV_corr_headers, on=(Observations.rootname \
                                             == NUV_corr_headers.rootname))
                    ).dicts())
    
    elif detector == 'FUV':
        data = list(
                    (Observations.select()
                    .join(FUVA_raw_headers, on=(Observations.rootname \
                                             == FUVA_raw_headers.rootname))
                    .join(FUVB_raw_headers, on=(Observations.rootname \
                                             == FUVB_raw_headers.rootname))
                    .join(FUVA_corr_headers, on=(Observations.rootname \
                                              == FUVA_corr_headers.rootname))
                    .join(FUVB_corr_headers, on=(Observations.rootname \
                                              == FUVB_corr_headers.rootname))
                    ).dicts())

    entry_to_print = filter(lambda dictionary: \
                            dictionary['rootname'] == rootname, data)
    
    database.close()
    print(entry_to_print)

def make_master_table():

    # Need to add:

    # PROP_TYP 
    # CAL_VER 
    # TAGFLASH (?)
    # FOCUS (Observations ?)
    # APER-DISP (Observations ?)
    # APER-XDISP (Observations ?)
    # LAPDSTP (Observations ?)
    # LAPXSTP (Observations ?)
    # OSM1-Coarse
    # OSM1-Fine
    # OSM2-Coarse
    # OSM2-Fine
    # DETHVL (NUV | rawtag)
    # DETHVLA (FUV | rawtag)
    # DETHVLB (FUV | rawtag)
    # OVERFLOW NUV CORR and FUV RAW
    
    settings = get_settings()
    
    database = get_database()
    database.connect()

    data = list(Observations.select(Observations.rootname, 
                                Observations.proposid,
                                Observations.exp_num, 
                                Observations.exptype,
                                Observations.targname,
                                Observations.extended,
                                Observations.obsmode,
                                Observations.postarg1,
                                Observations.postarg2,
                                Observations.obstype,
                                Observations.detector,
                                Observations.fppos,
                                Observations.cenwave,
                                Observations.aperture,
                                Observations.opt_elem,
                                Observations.date_obs,
                                Observations.expstart,
                                Observations.expend,
                                Observations.apmpos,
                                Observations.aperxpos,
                                Observations.aperypos,
                                Observations.opus_ver,
                                Observations.segment,
                                Observations.time_obs,
                                Observations.exptype,
                                Observations.asn_mtyp,
                                Observations.asn_id,
                                Observations.rawacq,
                                Observations.rawtag,
                                Observations.rawtag_a,
                                Observations.rawtag_b,
                                Observations.corrtag_a,
                                Observations.corrtag_b,
                                Observations.x1d_a,
                                Observations.x1d_b,
                                Observations.x1d,
                                Observations.csum_a,
                                Observations.csum_b,
                                Observations.flt, 
                                Observations.flt_a, 
                                Observations.flt_b, 
                                Observations.counts,
                                Observations.counts_a, 
                                Observations.counts_b, 
                                Observations.spt,
                                FUVA_raw_headers.hvlevela,
                                FUVA_raw_headers.neventsa,
                                FUVA_raw_headers.deventa,
                                FUVA_raw_headers.feventa,
                                FUVB_raw_headers.hvlevelb,
                                FUVB_raw_headers.neventsb,
                                FUVB_raw_headers.deventb,
                                FUVB_raw_headers.feventb).dicts()\
                            .join(FUVA_raw_headers, JOIN.LEFT_OUTER, on=(Observations.rootname == FUVA_raw_headers.rootname))\
                            .join(FUVB_raw_headers, JOIN.LEFT_OUTER, on=(Observations.rootname == FUVB_raw_headers.rootname))
        )

    database.close()
    df = pd.DataFrame(data)
    csv_path = os.path.join(settings['monitor_location'], 'master_table.csv')
    df.to_csv(csv_path)