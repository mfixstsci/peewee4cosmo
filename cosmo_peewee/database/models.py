from __future__ import print_function, absolute_import, division

import os

from peewee import *

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
    """ Parse config file and load settings

    If no config file is supplied, the configuration file will assume to be
    located at '~/configure.yaml'.

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
    database = MySQLDatabase(settings['database'],
                             host=settings['host'],
                             port=settings['port'],
                             user=settings['user'],
                             passwd=settings['password'])
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
    rootname = CharField()
    monitor_flag = BooleanField()

    class Meta:
        db_table = 'files'    
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

    filename = ForeignKeyField(Files,
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

    filename = ForeignKeyField(Files,
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

    filename = ForeignKeyField(Files,
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

    filename = ForeignKeyField(Files,
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

    filename = ForeignKeyField(Files,
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

    filename = ForeignKeyField(Files,
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

    filename = ForeignKeyField(Files,
                               db_column='filename',
                               default=None,
                               to_field="filename",
                               on_delete='CASCADE')
    class Meta:
        db_table = 'fuvb_corr_headers'

#-------------------------------------------------------------------------------
