"""Open pull and return header keywords for specific tables.

Some tables don't require any calculated calculated quantities. We just
pull and return the meta data that is populat4ed in the headers.
"""

import os

from astropy.io import fits
import gc

def nuv_corr_keys(file_result):
    """Keys for COS NUV corrtags
    
    Parameters
    ----------
    file_result: tuple
        A tuple with a path and filename
    
    Returns
    -------
    keywords: dict
        A dictionary that populates table metadata rows
    """

    file_path = os.path.join(file_result.path, file_result.filename)
    with fits.open(file_path) as hdu:
        keywords = {'filename': file_result.filename,
                    'rootname': hdu[0].header['rootname'],
                    'shift1a': hdu[1].header['shift1a'],
                    'shift1b': hdu[1].header['shift1b'],
                    'shift1c': hdu[1].header['shift1c'],
                    'shift1b': hdu[1].header['shift1b'],
                    'shift2b': hdu[1].header['shift2b'],
                    'shift2c': hdu[1].header['shift2c'],
                    'sp_loc_a': hdu[1].header['sp_loc_a'],
                    'sp_loc_b': hdu[1].header['sp_loc_b'],
                    'sp_loc_c': hdu[1].header['sp_loc_c'],
                    'sp_off_a': hdu[1].header['sp_off_a'],
                    'sp_off_b': hdu[1].header['sp_off_b'],
                    'sp_off_c': hdu[1].header['sp_off_c'],
                    'sp_nom_a': hdu[1].header['sp_nom_a'],
                    'sp_nom_b': hdu[1].header['sp_nom_b'],
                    'sp_nom_c': hdu[1].header['sp_nom_c'],
                    'sp_hgt_a': hdu[1].header['sp_hgt_a'],
                    'sp_hgt_b': hdu[1].header['sp_hgt_b'],
                    'sp_hgt_c': hdu[1].header['sp_hgt_c'],
                    'exptime': hdu[1].header['exptime']
                    }
    return keywords
    

def fuva_raw_keys(file_result):
    """Keys for COS segment FUVA rawtags
    
    Parameters
    ----------
    file_result: tuple
        A tuple with a path and filename
    
    Returns
    -------
    keywords: dict
        A dictionary that populates table metadata rows
    
    """

    file_path = os.path.join(file_result.path, file_result.filename)
    with fits.open(file_path) as hdu:
        keywords = {'filename': file_result.filename,
                    'rootname': hdu[0].header['rootname'],
                    'rawtime': hdu[1].header['rawtime'],
                    'neventsa': hdu[1].header['neventsa'],
                    'deventa': hdu[1].header['deventa'],
                    'feventa': hdu[1].header['feventa'],
                    'hvlevela': hdu[1].header['hvlevela']
                    }
    return keywords


def fuvb_raw_keys(file_result):
    """Keys for COS segment FUVB rawtags
    
    Parameters
    ----------
    file_result: tuple
        A tuple with a path and filename
    
    Returns
    -------
    keywords: dict
        A dictionary that populates table metadata rows
    
    """

    file_path = os.path.join(file_result.path, file_result.filename)
    with fits.open(file_path) as hdu:
        keywords = {'filename': file_result.filename,
                    'rootname': hdu[0].header['rootname'],
                    'rawtime': hdu[1].header['rawtime'],
                    'neventsb': hdu[1].header['neventsb'],
                    'deventb': hdu[1].header['deventb'],
                    'feventb': hdu[1].header['feventb'],
                    'hvlevelb': hdu[1].header['hvlevelb']
                    }
    return keywords


def fuva_corr_keys(file_result):
    """Keys for COS segment FUVA corrtags
    
    Parameters
    ----------
    file_result: tuple
        A tuple with a path and filename
    
    Returns
    -------
    keywords: dict
        A dictionary that populates table metadata rows
    
    """

    file_path = os.path.join(file_result.path, file_result.filename)
    with fits.open(file_path) as hdu:
        keywords = {'filename': file_result.filename,
                    'rootname': hdu[0].header['rootname'],
                    'shift1a': hdu[1].header['shift1a'],
                    'shift2a': hdu[1].header['shift2a'],
                    'sp_loc_a': hdu[1].header['sp_loc_a'],
                    'sp_off_a': hdu[1].header['sp_off_a'],
                    'sp_err_a': hdu[1].header.get('sp_err_a', -999.9),
                    'sp_nom_a': hdu[1].header['sp_nom_a'],
                    'sp_hgt_a': hdu[1].header['sp_hgt_a'],
                    'exptime': hdu[1].header['exptime']
                    }
    return keywords


def fuvb_corr_keys(file_result):
    """Keys for COS segment FUVB corrtags
    
    Parameters
    ----------
    file_result: tuple
        A tuple with a path and filename
    
    Returns
    -------
    keywords: dict
        A dictionary that populates table metadata rows
    
    """

    file_path = os.path.join(file_result.path, file_result.filename)
    with fits.open(file_path) as hdu:
        keywords = {'filename': file_result.filename,
                    'rootname': hdu[0].header['rootname'],
                    'shift1b': hdu[1].header['shift1b'],
                    'shift2b': hdu[1].header['shift2b'],
                    'sp_loc_b': hdu[1].header['sp_loc_b'],
                    'sp_off_b': hdu[1].header['sp_off_b'],
                    'sp_err_b': hdu[1].header.get('sp_err_b', -999.9),
                    'sp_nom_b': hdu[1].header['sp_nom_b'],
                    'sp_hgt_b': hdu[1].header['sp_hgt_b'],
                    'exptime': hdu[1].header['exptime']
                    }
    return keywords


def obs_keys(file_result):
    """Keys for observations
    
    Parameters
    ----------
    file_result: tuple
        A tuple with a path and filename
    
    Returns
    -------
    keywords: dict
        A dictionary that populates table metadata rows
    
    """

    file_path = os.path.join(file_result.path, file_result.filename)
    with fits.open(file_path) as hdu:
        keywords = {'rootname': hdu[0].header['rootname'],
                    'date_obs': hdu[1].header['date-obs'],
                    'detector': hdu[0].header['detector'],
                    'imagetyp': hdu[0].header['imagetyp'],
                    'targname': hdu[0].header['targname'],
                    'proposid': hdu[0].header['proposid'],
                    'ra_targ': hdu[0].header['ra_targ'],
                    'dec_targ': hdu[0].header['dec_targ'],
                    'pr_inv_l': hdu[0].header['pr_inv_l'],
                    'pr_inv_f': hdu[0].header['pr_inv_f'],
                    'opus_ver': hdu[0].header['opus_ver'],
                    'obstype': hdu[0].header['obstype'],
                    'obsmode': hdu[0].header['obsmode'],
                    'exptype': hdu[0].header['exptype'],
                    'postarg1': hdu[0].header['postarg1'],
                    'postarg2': hdu[0].header['postarg2'],
                    'life_adj': hdu[0].header['life_adj'],
                    'fppos': hdu[0].header['fppos'],
                    'exp_num': hdu[0].header['exp_num'],
                    'cenwave': hdu[0].header['cenwave'],
                    'propaper': hdu[0].header['propaper'],
                    'apmpos': hdu[0].header.get('apmpos', 'N/A'),
                    'aperxpos': hdu[0].header.get('aperxpos', -999.9),
                    'aperypos': hdu[0].header.get('aperypos', -999.9),
                    'aperture': hdu[0].header['aperture'],
                    'opt_elem': hdu[0].header['opt_elem'],
                    'extended': hdu[0].header['extended'],
                    'opus_ver': hdu[0].header['opus_ver'],
                    'life_adj': hdu[0].header['life_adj'],
                    'expstart': hdu[1].header['expstart'],
                    'expend': hdu[1].header['expend']
                    }
    return keywords


def file_keys(file_result):
    """Keys and tests for all data entering the system. This is so we can track 
    data files that have no data in them as well as avoid ingesting them into 
    any monitoring tables.
    
    Parameters
    ----------
    file_result: tuple
        A tuple with a path and filename
    
    Returns
    -------
    keywords: dict
        A dictionary that populates table metadata rows
    
    """
    path, filename = file_result
    
    try:
        # Set hard keys...
        keywords = {'path': path,
                    'filename': filename}

        # Open files to see if it opens or has data
        with fits.open(os.path.join(path,filename)) as hdu:
            if not len(hdu[1].data):
                keywords['monitor_flag'] = False
        
        # Hey, looks like it opens and has data
        keywords['monitor_flag'] = True
    
    # Sorry, file won't open or doesn't have data...
    except (IOError, TypeError, ValueError):
        keywords['monitor_flag'] = False
    
    # Add rootname to table
    rootname = filename.split('_')[0]
    if not len(rootname) == 9:
        rootname = None
    
    keywords['rootname'] = rootname

    return keywords 

def hv_keys(file_result):
    """Keys for the high voltage monitor
    
    Parameters
    ----------
    file_result: tuple
        A tuple with a path and filename
    
    Returns
    -------
    keywords: dict
        A dictionary that populates table metadata rows
    """

    file_path = os.path.join(file_result.path, file_result.filename)
    
    with fits.open(file_path) as hdu:
        segment = hdu[0].header['segment']
        keywords = {'filename':file_result.filename,
                    'rootname':hdu[0].header['rootname'],
                    'expstart':hdu[1].header['expstart'],
                    'segment':segment}
        
        if segment == 'FUVA':
            keywords['dethvl'] = hdu[1].header['dethvla']
            keywords['dethvc'] = hdu[1].header['dethvca']
            keywords['dethvn'] = hdu[1].header['dethvna']
            keywords['hvlevel'] = hdu[1].header['hvlevela']

        elif segment == 'FUVB':
            keywords['dethvl'] = hdu[1].header['dethvlb']
            keywords['dethvc'] = hdu[1].header['dethvcb']
            keywords['dethvn'] = hdu[1].header['dethvnb']
            keywords['hvlevel'] = hdu[1].header['hvlevelb']
    
    return keywords


def new_file_keys(file_result):
    path, filename = file_result

    try:
        # Set hard keys...
        keywords = {'path': path,
                    'filename': filename}

        # If no data, whoops, flag is false
        data = fits.getdata(os.path.join(path, filename))
        if not len(data):
            keywords['monitor_flag'] = False
        else:
            # Hey, looks like it opens and has data
            keywords['monitor_flag'] = True

        del data
        gc.collect()

    # Sorry, file won't open or doesn't have data...
    except (IOError, TypeError, ValueError, IndexError):
        keywords['monitor_flag'] = False

    # Add rootname to table
    rootname = filename.split('_')[0]
    if not len(rootname) == 9:
        rootname = None

    keywords['rootname'] = rootname

    return keywords
