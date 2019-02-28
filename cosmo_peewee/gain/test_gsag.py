import glob
from astropy.time import Time
import matplotlib
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 
import sys
sys.path.append('/user/duval/COS/common')
from cos_utils import *
from gain_monitoring import read_gsag_table
from gain_measurements import make_total_gainmap
import datetime

def test_gsag_tab_lp4_all_dates(segment, gsag_table_file = '/grp/crds/hst/references/hst/23e16470l_gsag.fits'):

    """ 
    This compares the GSAG table holes to the total gain maps every two weeks to track down when holes appear, since the beginning of LP4)
    """

    lp4_start = Time('2017-10-02', format = 'iso')
    lp4_start_mjd = lp4_start.mjd
    now = datetime.datetime.now()
    now_string = str(now.year) + '-'  + str(now.month) + '-' + str(now.day)
    nowt = Time(now_string, format = 'iso')
    now_mjd = nowt.mjd
    delta_t = 15. #interval between tests
    
    times = np.arange(lp4_start_mjd, now_mjd, delta_t)
    hvs = np.zeros_like(times)+ 163

    test_gsag_table(segment, this_time = times, this_hvs = hvs)
        


def test_gsag_table(segment, this_time=[], this_hvs=[], gsag_table_file = '/grp/crds/hst/references/hst/23e16470l_gsag.fits'):


    """
    This code tests that all gain sag holes are in teh gain sag table and that all gain sag holes in teh gain sag table are indeed gain sag holes in the gain maps. Writes a fits file with a binary mask 1/0 with 1  = hole, 0 = no hole, one mask from the total gain maps, one mask from teh gain sag table, and the difference between the two masks (should be 0).
    INPUTS:
    segment = "FUVA" or "FUVB"
    this time = a list of times to check gain sag for (will use a total gain map at those times)
    this_hvs  = a list of HV values for each of those times
    """

    

    if len(this_time)>0:
        time = this_time
        hvs = this_hvs
    else:
        delta_t = 30 #time resolution of the testing
        
        #create arrays of times and hv values over the history of COS/FUV
    
        tmin, tmax, hv = get_cos_hv_history(segment)
        tmin_mjd = Time(tmin, format = 'iso').mjd
        tmax_mjd = Time(tmax, format = 'iso').mjd
        
        time = []
        hvs = []

        for i in range(len(tmin)):
            t1 = tmin_mjd[i]
            t2 = tmax_mjd[i]
            this_hv = hv[i]
            this_time = list(np.arange(t1, t2, delta_t))
            this_hv= list(np.zeros(len(this_time)) + this_hv)
            time = time + this_time
            hvs = hvs + this_hv


        time = np.array(time)
        hvs = np.array(hvs)

    #this is the part of the detector we care about
    if segment=='FUVA':
        yrange=[400, 570]
    else:
        yrange=[460, 620]

    #for each time sample, we need to create a total gain map at each HV, and an array of 1/0 with 1 if it is a GSAG hole in the GSAG TAB. Then we compare the two.

    for (time, hv) in list(zip(time, hvs)):
        date = Time(time,format= 'mjd')
        date.format= 'datetime'
        year = date.value.year
        month = date.value.month
        day = date.value.day

        #read GSAG table 
        lx, ly, dx, dy, gsag_times = read_gsag_table(segment, hv, file=gsag_table_file)
        
        gain, counts, exptime=make_total_gainmap(hv, segment = segment, end_mjd = time, scale = True)

        #make a mask of the gain sag holes from the gain map
        mask_gain = np.zeros_like(gain)
        mask_gain[np.where(gain < 3.)] = 1.

        #make the mask from the GSAG table
        mask_gsag = make_gsag_mask(lx, ly, dx, dy, gsag_times, time)

        mask_dif  = mask_gain-mask_gsag

        outfile = '/grp/hst/cos/duval/TOTAL_GAIN/SCALED_TOTAL_GAIN_FOR_GSAG_TEST/gsag_test_{}_{}_{}-{}-{}'.format(segment, hv, year, month, day)

        hdu_out = fits.HDUList(fits.PrimaryHDU())

        # Adding primary header with file specifications to make results 
        # reproducible
        hdu_out[0].header['TELESCOP'] = 'HST'
        hdu_out[0].header['INSTRUME'] = 'COS'
        hdu_out[0].header['DETECTOR'] = 'FUV'
        hdu_out[0].header['OPT_ELEM'] = 'ANY'
        hdu_out[0].header['FILETYPE'] = 'GAINMAP'
        hdu_out[0].header['EXP_END'] = time 
        hdu_out[0].header['SEGMENT'] = segment
        hdu_out[0].header['HVLEVEL'] = hv

        # Data ext
       
        hdu_out.append(fits.ImageHDU(data=gain))
        hdu_out[1].header['EXTNAME'] = 'TOTAL GAIN'
        hdu_out.append(fits.ImageHDU(data=mask_gain))
        hdu_out[2].header['EXTNAME'] = 'MASK GAIN'
        hdu_out.append(fits.ImageHDU(data=mask_gsag))
        hdu_out[3].header['EXTNAME'] = 'MASK GSAG'
        hdu_out.append(fits.ImageHDU(data=mask_dif))
        hdu_out[4].header['EXTNAME'] = 'MASK DIF'

        
        hdu_out.writeto(outfile + ".fits", overwrite=True)
        hdu_out.close()
        

        

def make_gsag_mask(lx, ly, dx, dy, gsag_times, time):

    mask = np.zeros((1024, 16384), dtype = 'float32')

    holes_index = np.where(gsag_times <= time)
    holes_index = holes_index[0]

    glx = lx[holes_index]
    gly = ly[holes_index]
    gdx = dx[holes_index]
    gdy = dy[holes_index]

    for i in range(len(glx)):

        mask[gly[i]:gly[i] + gdy[i], glx[i]:glx[i] + gdx[i]] = 1.

    return(mask)
        
