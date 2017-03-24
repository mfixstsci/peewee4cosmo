from __future__ import print_function, absolute_import, division

import os

from peewee import *
from playhouse.pool import PooledMySQLDatabase

try:
    import yaml
except ImportError:
    from .yaml import yaml

__all__ = ['get_database']

#-------------------------------------------------------------------------------

def get_settings(config_file=None):
    """ Parse config file and load settings

    If no config file is supplied, the configuration file will assume to be
    located at '~/configure.yaml'.

    Parameters
    ----------
    config_file : str, optional
        yaml file containing configuration settings.

    Returns
    -------
    settings : dict
        dictionary of all settings

    """

    config_file = config_file or os.path.join(os.environ['HOME'], "configure_dev.yaml")
    
    with open(config_file, 'r') as f:
        settings = yaml.load(f)

    return settings
#-------------------------------------------------------------------------------

def get_database(config_file=None):
    """ Open connection to database

    Parameters
    ----------
    config_file : str, optional
        yaml file containing configuration settings.

    Returns
    -------
    Database : object
        A Peewee MySQLDatabase object.

    """

    settings = get_settings()

    database = PooledMySQLDatabase(settings['database'],
                                   host=settings['host'],
                                   port=settings['port'],
                                   user=settings['user'],
                                   passwd=settings['password'], stale_timeout=150)
    return database

#-------------------------------------------------------------------------------

class BaseModel(Model):
    """A base model that will use our MySQL database"""
    class Meta:
        database = get_database()

#-------------------------------------------------------------------------------

class Files(BaseModel):
    """Main files table"""
    path = CharField()
    filename = CharField(primary_key=True)
    monitor_flag = BooleanField()
    
    class Meta:
        db_table = 'files'   

#-------------------------------------------------------------------------------

class Observations(BaseModel):
    """Observations table"""
    path = CharField()
    filename = CharField(primary_key=True)
    targname = CharField()
    rootname = CharField()
    monitor_flag = BooleanField()
  
    class Meta:
        db_table = 'observations'  

#-------------------------------------------------------------------------------

class NUV_raw_headers(BaseModel):
    """NUV Raw Header Table"""

    rootname = CharField(primary_key=True)
    date_obs = CharField()
    detector = CharField()
    imagetyp = CharField()
    targname = CharField()
    proposid = FloatField()
    ra_targ = FloatField()
    dec_targ = FloatField()
    pr_inv_l = CharField()
    pr_inv_f = CharField()
    opus_ver = CharField()
    obstype = CharField()
    obsmode = CharField()
    exptype = CharField()
    postarg1 = IntegerField()
    postarg2 = IntegerField()
    life_adj = IntegerField()
    fppos = IntegerField()
    exp_num = IntegerField()
    cenwave = IntegerField()
    propaper = CharField()
    apmpos = CharField()
    aperxpos = FloatField()
    aperypos = FloatField()
    aperture = CharField()
    opt_elem = CharField()
    extended = CharField()
    obset_id = CharField()
    asn_id = CharField()
    asn_mtyp = CharField()
    expstart = FloatField()
    expend = FloatField()
    exptime = FloatField()
    nevents = IntegerField()

    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')

    class Meta:
        db_table = 'nuv_raw_headers'
 
#-------------------------------------------------------------------------------

class NUV_corr_headers(BaseModel):
    """NUV Corrtag Header Table"""

    rootname = CharField(primary_key = True)
    shift1a = FloatField()
    shift1b = FloatField()
    shift1c = FloatField()
    shift1b = FloatField()
    shift2b = FloatField()
    shift2c = FloatField()
    sp_loc_a = FloatField()
    sp_loc_b = FloatField()
    sp_loc_c = FloatField()
    sp_off_a = FloatField()
    sp_off_b = FloatField()
    sp_off_c = FloatField()
    sp_nom_a = FloatField()
    sp_nom_b = FloatField()
    sp_nom_c = FloatField()
    sp_hgt_a = FloatField()
    sp_hgt_b = FloatField()
    sp_hgt_c = FloatField()
    exptime = FloatField()

    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')

    class Meta:
        db_table = 'nuv_corr_headers'
#-------------------------------------------------------------------------------

class FUV_primary_headers(BaseModel):
    """Table of shared keywords for both segments"""

    rootname = CharField(primary_key = True)
    date_obs = CharField()
    detector = CharField()
    imagetyp = CharField()
    targname = CharField()
    proposid = FloatField()
    ra_targ = FloatField()
    dec_targ = FloatField()
    pr_inv_l = CharField()
    pr_inv_f = CharField()
    opus_ver = CharField()
    obstype = CharField()
    obsmode = CharField()
    exptype = CharField()
    postarg1 = IntegerField()
    postarg2 = IntegerField()
    life_adj = IntegerField()
    fppos = IntegerField()
    exp_num = IntegerField()
    cenwave = IntegerField()
    propaper = CharField()
    apmpos = CharField()
    aperxpos = FloatField()
    aperypos = FloatField()
    aperture = CharField()
    opt_elem = CharField()
    extended = CharField()
    obset_id = CharField()
    asn_id = CharField()
    asn_mtyp = CharField()

    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    class Meta:
        db_table = 'fuv_primary_headers'

#-------------------------------------------------------------------------------

class FUVA_raw_headers(BaseModel):
    'Table of raw fuva data'

    rootname = CharField(primary_key = True)
    expstart = FloatField()
    expend = FloatField()
    rawtime = FloatField()
    neventsa = FloatField()
    deventa = FloatField()
    feventa = FloatField()
    hvlevela = IntegerField()

    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    class Meta:
        db_table = 'fuva_raw_headers'

#-------------------------------------------------------------------------------

class FUVB_raw_headers(BaseModel):
    'Table of raw fuvb data'

    rootname = CharField(primary_key = True)
    expstart = FloatField()
    expend = FloatField()
    rawtime = FloatField()
    neventsb = FloatField()
    deventb = FloatField()
    feventb = FloatField()
    hvlevelb = IntegerField()

    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    class Meta:
        db_table = 'fuvb_raw_headers'

#-------------------------------------------------------------------------------

class FUVA_corr_headers(BaseModel):
    'Table of corr fuva data'

    rootname = CharField(primary_key = True)
    shift1a = FloatField()
    shift2a = FloatField()
    sp_loc_a = FloatField()
    sp_off_a = FloatField()
    sp_err_a = FloatField()
    sp_nom_a = FloatField()
    sp_hgt_a = IntegerField()
    exptime = FloatField()

    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    class Meta:
        db_table = 'fuva_corr_headers'

#-------------------------------------------------------------------------------

class FUVB_corr_headers(BaseModel):
    'Table of corr fuvb data'

    rootname = CharField(primary_key = True)
    shift1b = FloatField()
    shift2b = FloatField()
    sp_loc_b = FloatField()
    sp_off_b = FloatField()
    sp_err_b = FloatField()
    sp_nom_b = FloatField()
    sp_hgt_b = IntegerField()
    exptime = FloatField()

    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    class Meta:
        db_table = 'fuvb_corr_headers'

#-------------------------------------------------------------------------------

class Lampflash(BaseModel):

    """Preped metadata for OSM monitor."""

    rootname = CharField(default='N/A')
    date = FloatField(default=-999.9)
    proposid = IntegerField(default=-999)
    detector = CharField(default='N/A')
    segment = CharField(default='N/A')
    opt_elem = CharField(default='N/A')
    cenwave = IntegerField(default=-999)
    fppos = IntegerField(default=-999)
    lamptab = CharField(default='N/A')
    flash = IntegerField(default=-999)
    x_shift = FloatField(default=-999.9)
    y_shift = FloatField(default=-999.9)
    filetype = CharField(default='N/A')
    cal_date = CharField(default='N/A')
    found = BooleanField(default=False)
    
    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    
    class Meta:
        db_table = 'lampflash'
#-------------------------------------------------------------------------------

class Rawacqs(BaseModel):

    """Preped metadata for OSM monitor."""
    
    #rootname = CharField(primary_key=True)
    rootname = CharField()
    date = FloatField()
    proposid = IntegerField()
    detector = CharField()
    opt_elem = CharField()
    cenwave = IntegerField()
    fppos = IntegerField()
    flash = IntegerField()
    x_shift = FloatField()
    y_shift = FloatField()
    filetype = CharField()
    cal_date = CharField()
    found = BooleanField()
    
    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    
    class Meta:
        db_table = 'rawacqs'

#-------------------------------------------------------------------------------

class Darks(BaseModel):
    """Record dark rate"""
    
    rootname = CharField()
    detector = CharField()
    date = FloatField()
    dark = FloatField()
    ta_dark = FloatField()
    latitude = FloatField()
    longitude = FloatField()
    sun_lat = FloatField()
    sun_lon = FloatField()
    temp = FloatField()

    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    
    class Meta:
        db_table = 'darks'

#-------------------------------------------------------------------------------

class Stims(BaseModel):
    """Record location of all STIM pulses"""

    rootname = CharField()
    time = FloatField(default=0)
    abs_time = FloatField(default=0)
    stim1_x = FloatField(default=0)
    stim1_y = FloatField(default=0)
    stim2_x = FloatField(default=0)
    stim2_y = FloatField(default=0)
    counts = FloatField(default=0)
    segment = CharField()
  
    filename = ForeignKeyField(Observations,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')

    class Meta:
        db_table = 'stims'
    
#-------------------------------------------------------------------------------