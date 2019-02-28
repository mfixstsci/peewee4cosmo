import sys
sys.path.append("PEEWEE4COSMO/cosmo_peewee/database")
sys.path.append("/user/duval/COS/common")
from cos_utils import mk_cos_image_jrd
from models import get_database, Observations, FUVB_raw_headers, FUVA_raw_headers
from astropy.time import Time
from astropy.table import Table
import datetime
from astropy.io import fits, ascii
from astropy.modeling import models, fitting
import numpy as np
from astropy.modeling import models, fitting
import os
import math
import time
import glob
from gain_monitoring import resample_gain, get_cos_hv_history
from matplotlib import pyplot as plt
import subprocess



def convert_figures():

    #files = glob.glob('/grp/hst/cos/duval/TOTAL_GAIN/SCALED_TOTAL_GAIN_FOR_MOVIE/map_gain*.png')
    files = glob.glob('/grp/hst/cos/duval/TOTAL_GAIN/SCALED_TOTAL_GAIN_FOR_MOVIE/IMG-*.png')
    
    for file in files:

        dir, root =os.path.split(file)
        dum = root.split('-')
        index = int(dum[len(dum)-1].replace('.png', ''))
        #print(index)
        new_file = dir + '/' + 'IMG-' + '{0:04d}'.format(index) + '.png'
        print(file, ' ', new_file)

        #new_file = file.replace('0', '')
        
        subprocess.call(['cp', file, new_file])


def plot_pixel_history(x, y, lp, segment):

    c = fits.open("LP{}_{}_cube_gain_counts.fits".format(lp, segment))
    gain = c[1].data
    counts = c[2].data
    time = c[3].data
    hv = c[4].data

    s = gain.shape
    nt= s[0]
    ny = s[1]
    nx = s[2]

    delta_counts = np.zeros((nt-1, ny, nx), dtype  = 'float32')
    delta_gain = np.zeros((nt-1, ny, nx), dtype  = 'float32')

    y0 = c[0].header['Y0']


    for i in range(1,nt-1,1):
        delta_counts[i-1, :, :] = counts[i,:,:]  - counts[i-1,:,:]
        delta_gain[i-1,:,:] = gain[i,:,:]  - gain[i-1,:,:]


    plt.plot(time[1:], delta_counts[:,y-y0,x], 'ko')
    plt.plot(time[1:], delta_gain[:, y-y0, x], 'ro')
    plt.show()
    
    return(time, y0, delta_counts, delta_gain)
        
        
    
def make_movie_gain(segment, delta_t = 15., redo=False):


    """
    This function tracks down the HV history of COS and creates maps of total gain adn total counts every delta_t for the entire history of COS
    """

    tmin, tmax, hv = get_cos_hv_history(segment)
    tmin_mjd = Time(tmin, format = 'iso').mjd
    tmax_mjd = Time(tmax, format = 'iso').mjd

    #since tmin starts at the begninig of COS, I want the fiest time to be at tmin + delta_t since it will look at gain maps before that

    tmin_mjd[0] = tmin_mjd[0] + delta_t

    print("HISTORY")

    print(tmin)
    print(tmax)
    print(tmin_mjd)
    print(tmax_mjd)
    print(hv)

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


    time = np.array(time, dtype = 'float32') 
    hvs = np.array(hvs, dtype = 'uint32')


    nt = len(time)

    #total_counts = np.zeros((nt, 1024, 16384), dtype = 'float32')
    #total_exptime = np.zeros(nt, dtype = 'float32')

    

    counter=0

    for (time, hv) in list(zip(time, hvs)):
        print(time, hv)
        date = Time(time,format= 'mjd')
        date.format= 'datetime'
        year = date.value.year
        month = date.value.month
        day = date.value.day

        fitsname = '/grp/hst/cos/duval/TOTAL_GAIN/SCALED_TOTAL_GAIN_FOR_MOVIE/total_gain_{}_{}_{}-{}-{}'.format(segment, hv, year, month, day)
        pngname = "/grp/hst/cos/duval/TOTAL_GAIN/SCALED_TOTAL_GAIN_FOR_MOVIE/IMG-{}-{}.png".format(segment, counter)

        if (os.path.isfile(fitsname)==False) or (redo==True):
        
            gain,counts, exptime=make_total_gainmap(hv, segment = segment, end_mjd = time, scale = True, writefits = fitsname, writepng = pngname)

            #total_counts[counter, :, :] = counts #total_counts[counter-1, :, :] + counts
            #total_exptime[counter] = exptime
            
            #if counter >=1:
            #I am making my own, but instead of going back to corrtags, I could go through all CCI in teh given time interval and add the counts. Actually going through the corrtafs does not work at LP4 due to proprietary data! Need to change hte total gain map code to count the counts.
            #counts, exp_time,dum_gain, pha_axis = create_cci(1234, hv, 'all', 1234, segment,time-delta_t ,time, compute_gain=False, cube=False)
            #counts = counts[0,:,:]
            

            #hdu  = fits.open('/grp/hst/cos/duval/TOTAL_GAIN/SCALED_TOTAL_GAIN_FOR_MOVIE/total_gain_{}_{}_{}-{}-{}_scaled'.format(segment, hv, year, month, day) + '.fits', mode = 'update')
            #nhdu = len(hdu)
            #hdu.append(fits.ImageHDU(data=total_counts[counter, :, :]))
            #hdu[nhdu].header['EXTNAME'] = 'TOTAL COUNTS'

            #hdu.flush()
            #hdu.close()

            
        counter+=1
        





def make_total_gainmap(hv_lvl, **kwargs):
    """Make total gainmaps for each HV level and the total over gainmap.
    Parameters
    ----------
    hv_lvl: int
        High voltage for gainmap you want to make.
    **kwargs
        Arbitrary number of keyword arguements.
    
    Returns
    -------
    gainmap
    """

    
    gainmap_dir = kwargs.get('gainmap_dir', '/grp/hst/cos/Monitors/CCI/')
    segment = kwargs.get('segment', 'FUVB')
    start_mjd = kwargs.get('start_mjd', 55055)
    end_mjd = kwargs.get('end_mjd', 70000)
    reverse = kwargs.get('reverse', False)
    writefits = kwargs.get('writefits', '')
    writepng = kwargs.get('writepng', '')
    resamp = kwargs.get('resamp', True)
    mincounts = kwargs.get('mincounts', 50)
    scale= kwargs.get('scale', False)
    cmap = kwargs.get('cmap', 'magma')
    save_pix_hist = kwargs.get('save_pix_hist', False)
    x =  kwargs.get('x', 0)
    y = kwargs.get('y', 0)

    count_hist = []
    gain_hist = []
    time_hist = []

    # Depending on the segment, the filename are different.

    if scale==False:
        if segment == 'FUVA':
            search_string = 'l_*_00_{}_cci_gainmap.fits*'.format(hv_lvl)
        elif segment == 'FUVB':
            search_string = 'l_*_01_{}_cci_gainmap.fits*'.format(hv_lvl)
    else:
        if segment == 'FUVA':
            search_string = 'l_*_00_*_cci_gainmap.fits*'
        elif segment == 'FUVB':
            search_string = 'l_*_01_*_cci_gainmap.fits*'
            
    # Get all of the data and sort.
    all_datasets = [item for item in glob.glob(os.path.join(gainmap_dir, 
                                                            search_string))]
    all_datasets.sort()
    # Switch the order to reverse to see what the detector originally looked 
    # like.
    if reverse:
        all_datasets = all_datasets[::-1]

    # Make an empty array the size of the detector.
    out_gain = np.zeros((512,2048)) + np.nan
    out_counts= np.zeros((512,2048))
    out_exptime = 0.


    # For each gainmap
    for item in all_datasets:
        
        cci_hdu = fits.open(item)

    
        # Make sure you are making a total gainmap in the dates you defined.
        if not cci_hdu[0].header['EXPSTART'] >= start_mjd:
            continue
        
        if not cci_hdu[0].header['EXPSTART'] <= end_mjd:
            continue     
        
        # Get the data
        det_hv = cci_hdu[0].header['DETHV']
        cci_gain_data = cci_hdu['MOD_GAIN'].data
        cci_count_data = cci_hdu['COUNTS'].data

        if det_hv != hv_lvl:
            delta_gain = 0.393*(det_hv-hv_lvl)
            cci_gain_data = cci_gain_data-delta_gain
            sag = np.where(cci_gain_data <0.)
            if len(sag[0]) > 0:
                cci_gain_data[sag] = 0.
            #print("CHECK HV ", hv_lvl, det_hv, delta_gain)

        #print("DATE ")
        #print(item, cci_hdu[0].header['EXPSTART'])
        #print("GAIN ", cci_gain_data[290, 1108])
        #print("COUNTS ", cci_count_data[290, 1108])

        #JRD: This does not work. In the middle of gain sag holes, the gain is measured as 0.0, the same as pixels with no coverage. As a result, the gain in the center of the gain sag holes is not updated, leading to high gain values there. 
        # Find where there is data 
        #index = np.where(cci_data)

        #JRD FIx to problem above
        index = np.where(cci_count_data >=mincounts )
        
        # get the high voltage level.
        # dethv = cci_hdu[0].header['DETHV']

        # Add the data to the array and move to the next gainmap.
        out_gain[index] = cci_gain_data[index]
        out_counts[index] = out_counts[index] + cci_count_data[index]
        out_exptime += cci_hdu[0].header['EXPTIME']

        if save_pix_hist==True:

            #if cci_count_data[495//2, 5052//8] >= mincounts:
            #    print(item)
            #    print(cci_count_data[495//2, 5052//8])

            if resamp==True:
                resamp_cci_gain = resample_gain('', 'MOD_GAIN', data = cci_gain_data)
                resamp_cci_counts = resample_gain('', 'COUNTS', data = cci_count_data)/16. #converting counts in an 8x2 superpixel to 1 pixel
            else:
                resamp_cci_gain = cci_gain_data
                resamp_cci_counts = cci_count_data


            #Save the history of a pixel
            if resamp_cci_counts[y,x] >= mincounts:
                count_hist.append(resamp_cci_counts[y,x])
                gain_hist.append(resamp_cci_gain[y,x])
                time_hist.append(cci_hdu[0].header['EXPSTART'])
                              

    if resamp==True:
        resamp_gain = resample_gain('', 'MOD_GAIN', data = out_gain)
        resamp_counts = resample_gain('', 'COUNTS', data = out_counts)/16.
        
    else:
        resamp_gain = out_gain
        resamp_counts = out_counts
        

                #save the history a pixel

    if writefits != '':

        #hdu = fits.HDUList(fits.PrimaryHDU())
        #hdu.append(fits.ImageHDU(data = out_counts))
        #hdu[1].header['EXTNAME']= 'COUNTS ORIG SIZE'
        #hdu.writeto('test_counts.fits', overwrite=True)

        
        filename = writefits.replace(".fits", "")
        if scale==True:
            filename = filename + "_scaled"
        hdu_out = fits.HDUList(fits.PrimaryHDU())

        # Adding primary header with file specifications to make results 
        # reproducible
        hdu_out[0].header['TELESCOP'] = 'HST'
        hdu_out[0].header['INSTRUME'] = 'COS'
        hdu_out[0].header['DETECTOR'] = 'FUV'
        hdu_out[0].header['OPT_ELEM'] = 'ANY'
        hdu_out[0].header['FILETYPE'] = 'GAINMAP'
                
        hdu_out[0].header['EXPSTART'] = start_mjd
        hdu_out[0].header['EXP_END'] = end_mjd
        hdu_out[0].header['CCI_DIR'] = gainmap_dir
        hdu_out[0].header['SEGMENT'] = segment
        hdu_out[0].header['HVLEVEL'] = hv_lvl
        hdu_out[0].header['DATEMADE'] = (time.strftime("%d/%m/%Y"))
        hdu_out[0].header['EXPTIME'] = out_exptime

        # Data ext
        
        hdu_out.append(fits.ImageHDU(data=resamp_gain))
        hdu_out[1].header['EXTNAME'] = 'TOTAL GAIN'

        hdu_out.append(fits.ImageHDU(data=resamp_counts))
        hdu_out[2].header['EXTNAME'] = 'TOTAL COUNTS'
        
        
        hdu_out.writeto(filename + ".fits", overwrite=True)
        hdu_out.close()

    if writepng != '':

        if segment=='FUVA':
            yrange=[400, 570]
        else:
            yrange=[460, 620]

        date = Time(end_mjd,format= 'mjd')
        date.format= 'datetime'
        year = date.value.year
        month = date.value.month
        day = date.value.day
        
        plt.clf()
        plt.close()
        fig,ax = plt.subplots(2,sharex = True,figsize=(15, 15))
        cax0=ax[0].imshow(resamp_gain, vmin = 0, vmax = 12, cmap = cmap, aspect = 'auto', origin = 'lower')
        ax[0].set_ylim(bottom = yrange[0], top = yrange[1])
        ax[0].text(8000, 600, '{} {} {}-{}-{}'.format(segment, hv_lvl, year, month, day), fontsize = 15, color = 'magenta')
        ax[0].set_xlabel("XCORR", fontsize = 18)
        ax[0].set_ylabel("YCORR", fontsize = 18)
        cbar0 = fig.colorbar(cax0, ax = ax[0])
        cbar0.ax.set_ylabel("MODAL GAIN", fontsize = 18)

        cax1=ax[1].imshow(resamp_counts, vmin = 0, vmax = 3.e4, cmap = 'magma', aspect = 'auto', origin = 'lower')
        ax[1].set_ylim(bottom = yrange[0], top = yrange[1])
        ax[1].text(8000, 600, '{} {} {}-{}-{}'.format(segment, hv_lvl, year, month, day), fontsize = 15, color = 'magenta')
        ax[1].set_xlabel("XCORR", fontsize = 18)
        ax[1].set_ylabel("YCORR", fontsize = 18)
        cbar1 = fig.colorbar(cax1, ax = ax[1])
        cbar1.ax.set_ylabel("TOTAL COUNTS", fontsize = 18)

        fig.subplots_adjust(hspace = 0.15, wspace =0, top = 0.95, bottom =0.12, left = 0.12, right = 0.95)

        plt.savefig(writepng , format = 'png', dpi = 1000)
        plt.clf()
        plt.close()
        

    if save_pix_hist==True:
        return(resamp_gain, resamp_counts, out_exptime, count_hist, gain_hist, time_hist)
    else:
        return(resamp_gain, resamp_counts, out_exptime)


def make_some_gainmaps(times, hvs, segments):

    for (time, hv) in list(zip(times, hvs)):
        for segment in segments:

            if type(time)==str:
                date = Time(time, format = 'iso')
                time_mjd = date.mjd
                date.format = 'datetime'
                date = str(date.value.year) + '-' + str(date.value.month) + '-' + str(date.value.day)
                
            else:    
                date = Time(time, format = 'mjd')
                date.format= 'datetime'
                date = str(date.value.year) + '-' + str(date.value.month) + '-' + str(date.value.day)
                time_mjd = time
            fitsname = 'total_gain_{}_{}_{}'.format(segment, hv, date)
            pngname = fitsname + '.png'
            gain,counts, exptime=make_total_gainmap(hv, segment = segment, end_mjd = time_mjd, scale = True, writefits = fitsname, writepng = pngname)
        
    

def make_all_gainmaps():
    """Make all of the total gainmaps.
    Parameters
    ----------
    hv_lvl : integer
        Command high voltage level
    **kwargs
        An arbitrary number of keyword arguements
    
    Returns
    -------
    None
    """
    gainmap_dir = '/grp/hst/cos/Monitors/CCI/'
    outdir = '/grp/hst/cos/duval/TOTAL_GAIN/'
    
    start_mjd =  Time('2009-08-12').mjd
    delta_t = 30
    now = Time(datetime.datetime.now()).mjd
    end_mjds =   np.arange(start_mjd +delta_t, now, delta_t)

    hv_lvls_ab = [[163, 167, 169, 171, 175, 178],[163, 167, 169, 175]]
    segments = ['FUVA', 'FUVB']

    for iseg in range(2):
        hv_lvls=hv_lvls_ab[iseg]
        segment = segments[iseg]

        for hv_lvl in hv_lvls:
        
            for end_mjd in end_mjds:

                this_end_mjd =Time(end_mjd, format = 'mjd')
                end_full = Time(this_end_mjd, format = 'iso', out_subfmt = 'date')
                end = end_full.value
                end_fits = Time(this_end_mjd, format = 'fits')
            
                filename = os.path.join(outdir,'total_gain_{}_{}_{}.fits'.format(segment,hv_lvl, end))
                hdu_out = fits.HDUList(fits.PrimaryHDU())

                # Adding primary header with file specifications to make results 
                # reproducible
                hdu_out[0].header['TELESCOP'] = 'HST'
                hdu_out[0].header['INSTRUME'] = 'COS'
                hdu_out[0].header['DETECTOR'] = 'FUV'
                hdu_out[0].header['OPT_ELEM'] = 'ANY'
                hdu_out[0].header['FILETYPE'] = 'GAINMAP'
                
                hdu_out[0].header['EXPSTART'] = start_mjd
                hdu_out[0].header['EXP_END'] = end_mjd
                hdu_out[0].header['CCI_DIR'] = gainmap_dir
                hdu_out[0].header['SEGMENT'] = segment
                hdu_out[0].header['HVLEVEL'] = hv_lvl
                hdu_out[0].header['DATEMADE'] = (time.strftime("%d/%m/%Y"))

                # Data ext
       
                hdu_out.append(fits.ImageHDU(data=make_total_gainmap(hv_lvl, 
                                                             gainmap_dir=gainmap_dir, 
                                                             segment=segment, 
                                                             start_mjd=start_mjd, 
                                                             end_mjd=end_mjd)))
                hdu_out[1].header['EXTNAME'] = 'TOTAL GAIN'
        
        
        
                hdu_out.writeto(filename, overwrite=True)
                hdu_out.close()


def play_w_db():

    all_obs = Observations.select().limit(100)
    fuvb4 = Observations.select().where((Observations.life_adj==4) & (Observations.cenwave==1291))

    db = get_database()
    db.connect()

    for obs in fuvb4:
        print(obs.life_adj)
        print(obs.detector)
        #print(obs.corrtag_b)
        print(obs.segment)
        print(obs.rootname)
        print(obs.proposid)
        print(obs.time_obs)
        print(obs.exptime)
        print(obs.expstart)
        print(obs.date_obs)

    db.close()

def create_cci(life_adj, hv, cenwave, fppos, segment, start_mjd, stop_mjd, compute_gain=False, cube=True,writefits = '', pha_low = 1, pha_high = 31, list_obs = False):

    if segment =='FUVB':
        ss = 'b'
    else:
        ss = 'a'

    ny = 1024
    nx = 16384
    
    #Create the arrays that will hold the counts, gain and exposure time

    #pha_low = 1
    #pha_high = 31
    npha = pha_high-pha_low + 2
    pha_axis = np.arange(pha_low, pha_high + 1, 1)

    if cube==True:
        counts = np.zeros((npha, ny, nx), dtype = 'float32')
    else:
        counts = np.zeros((ny,nx), dtype = 'float32')
        #if compute_gain==True:
    gain = np.zeros((ny, nx),dtype = 'float32')
    exp_time = 0.

    pids = []
    rootnames = []
    asns =[]


    #With all HV settings:

    if hv=='all':

        if life_adj == 1234:
            if segment =='FUVB':
                fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV')  & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
            else:
                fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))

        else:
            if cenwave=='all':
                if segment =='FUVB':
                    fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (Observations.life_adj==life_adj) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
                else:
                    fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.life_adj==life_adj) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))

            else:
    
                if fppos ==1234:
                    if segment =='FUVB':
                        fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.life_adj==life_adj) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
                    else:
                        fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.life_adj==life_adj) & (Observations.expstart >=start_mjd) & (Observations.expstart < stop_mjd))
       
                else:
                    if segment =='FUVB':
                        fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.fppos ==fppos) & (Observations.life_adj==life_adj) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
                    else:
                        fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.fppos ==fppos) & (Observations.life_adj==life_adj) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))


    #With a given HV setting

    else:

        if life_adj == 1234:
            if segment =='FUVB':
                fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (FUVB_raw_headers.hvlevelb==hv) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
            else:
                fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (FUVA_raw_headers.hvlevela==hv)& (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))

        else:
            if cenwave=='all':
                if segment =='FUVB':
                    fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (Observations.life_adj==life_adj) & (FUVB_raw_headers.hvlevelb==hv) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
                else:
                    fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.life_adj==life_adj) & (FUVA_raw_headers.hvlevela==hv)& (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))

            else:
    
                if fppos ==1234:
                    if segment =='FUVB':
                        fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.life_adj==life_adj) & (FUVB_raw_headers.hvlevelb==hv) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
                    else:
                        fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.life_adj==life_adj) & (FUVA_raw_headers.hvlevela==hv)& (Observations.expstart >=start_mjd) & (Observations.expstart < stop_mjd))
       
                else:
                    if segment =='FUVB':
                        fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.fppos ==fppos) & (Observations.life_adj==life_adj) & (FUVB_raw_headers.hvlevelb==hv) & (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
                    else:
                        fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.fppos ==fppos) & (Observations.life_adj==life_adj) & (FUVA_raw_headers.hvlevela==hv)& (Observations.expstart >= start_mjd) & (Observations.expstart < stop_mjd))
    
        

    print("FOUND {} {}  Observations taken between {} and {} with cenwave {} FPPOS {} at lifetime postion {} and HV-{} = {}".format(len(fuv_obs), segment,start_mjd, stop_mjd, cenwave, fppos, life_adj, segment, hv))

    if len(fuv_obs)>0:

        for obs in fuv_obs:
            rootname = obs.rootname
            asn = obs.asn_id
            pid = int(obs.proposid)
            path = '/grp/hst/cos2/cosmo/{}/{}_corrtag_{}.fits.gz'.format(np.round(pid), rootname, ss)
            print(pid, rootname)
            pids.append(pid)
            rootnames.append(rootnames)
            asns.append(asn)
            
            if list_obs ==False:
                
                #For testing purposes, I am going to ignore proprietary data
                if os.access(path, os.R_OK):
                    corrtag = fits.open(path)
                    ctag = corrtag[1].data
                    if cube==True:
                        cube = make_cube_pha(ctag, phamin = pha_low, phamax = pha_high)
                        counts[:,:,:] = counts[:,:,:] + cube
                    else:
                        these_counts, bpix, gsag_im = mk_cos_image_jrd('', corrtag_in = ctag,  pha_filter = False, no_flat =True)
                        counts[:,:]= counts[:,:] + these_counts
                        exp_time = exp_time +  obs.exptime
                else:
                    print("YOU DO NOT HAVE PERMISSION TO OPEN FILE ", path)
                if compute_gain==True:
                    gain[:,:], gstd = measure_gainimage(counts[1:,:,:], phlow = pha_low, phhigh = pha_high)
            else:
                print(pid, rootname, asn)

    if (writefits != '') and (list_obs==False):

        filename = writefits.replace('.fits', '')
        if (cenwave==1280) or (cenwave==800) or (cenwave==1105) or (cenwave==1105):
            opt_elem = 'G140L'
        if (cenwave ==1055) or (cenwave==1096) or (cenwave==1222) or (cenwave==1291) or (cenwave==1300) or (cenwave==1309) or (cenwave ==1318) or (cenwave==1327):
            opt_elem = 'G130M'
        if (cenwave ==1533) or (cenwave==1577) or (cenwave==1589) or (cenwave==1600) or (cenwave==1611) or (cenwave==1623):
            opt_elem = 'G160M'
        if cenwave == 'all':
            opt_elem = 'ALL'

        hduout = fits.HDUList(fits.PrimaryHDU())
        hduout[0].header['TELESCOP'] = 'HST'
        hduout[0].header['DETECTOR'] = 'FUV'
        hduout[0].header['OPT_ELEM'] = opt_elem
        hduout[0].header['CENWAVE'] = cenwave
        hduout[0].header['FPPOS']=fppos
        hduout[0].header['SEGMENT'] = segment
        hduout[0].header['HVLEVEL'] = hv
        hduout[0].header['FILETYPE'] = 'JRDCCI'
        hduout[0].header['EXPSTART'] = start_mjd
        hduout[0].header['EXP_END'] = stop_mjd
        hduout[0].header['EXPTIME'] = exp_time
        
        
        hduout.append(fits.ImageHDU(counts))
        hduout[1].header['EXTNAME'] = 'COUNTS'

        if cube==True:
            hduout.append(fits.ImageHDU(pha_axis))
            hduout[2].header['EXTNAME'] = 'PHA'
        hduout.writeto(filename+ '_CCI.fits', overwrite=True)
        
        if compute_gain==True:
            hdug = fits.HDUList(fits.PrimaryHDU())
            hdug[0].header['TELESCOP'] = 'HST'
            hdug[0].header['DETECTOR'] = 'FUV'
            hdug[0].header['OPT_ELEM'] = opt_elem
            hdug[0].header['CENWAVE'] = cenwave
            hdug[0].header['FPPOS']=fppos
            hdug[0].header['SEGMENT'] = segment
            hdug[0].header['HVLEVEL'] = hv
            hdug[0].header['FILETYPE'] = 'JRDCCI'
            hdug[0].header['EXPSTART'] = start_mjd
            hdug[0].header['EXP_END'] = stop_mjd
            hdug.append(fits.ImageHDU(gain))
            hdug[1].header['EXTNAME'] = 'MOD_GAIN'
            hdug.writeto(filename+ '_GAIN.fits', overwrite=True)
            
    if list_obs==False:
        return(counts, exp_time,gain, pha_axis)
    else:
        return(pids, rootnames, asns)


def make_cube_pha(corrtag, phamin = 1, phamax = 31):

    xmin = 0
    xmax = 16383
    ymin=0
    ymax=1023
    ny = 1024
    nx = 16384
    npha = phamax - phamin  +2  #for total counts
    cube = np.zeros((npha, ny, nx), dtype = 'float32')

    xcorr = corrtag['XCORR']
    ycorr = corrtag['YCORR']
    pha  =corrtag['PHA']

    valid= np.where((xcorr >= 0) & (xcorr < 16383) & (ycorr >=0) & (ycorr < 1023))
    valid = valid[0]
    
    xbin = np.asarray(np.floor((xcorr[valid] + 0.5)), dtype=np.int)
    ybin = np.asarray(np.floor((ycorr[valid] + 0.5)), dtype=np.int)
    phabin=np.asarray(np.floor((pha[valid] + 0.5)), dtype=np.int)

    for i in range(len(xbin)):
        if (phabin[i] >= phamin) and (phabin[i] <= phamax):
            cube[phabin[i]-phamin+1, ybin[i], xbin[i]]+=1.

    cube[0,:,:] = np.sum(cube[1:,:,:], axis = 0)

    return(cube)


def measure_gainimage(data_cube, **kwargs):
    """Measure the modal gain at each pixel
    returns a 2d gainmap.
    """

    mincounts = kwargs.get('mincounts',10)
    phlow = kwargs.get('phlow', 1)
    phhigh = kwargs.get('phhigh', 31)
    pha_axis = np.arange(phlow, phhigh+1, 1)
    plot = kwargs.get('plot', False)

    # Suppress certain pharanges
    #for i in list(range(0, phlow+1)) + list(range(phhigh, len(data_cube))):
    #    data_cube[i] = 0

    counts_im = np.nansum(data_cube, axis=0)

    out_gain = np.zeros(counts_im.shape)
    out_counts = np.zeros(counts_im.shape)
    out_std = np.zeros(counts_im.shape)

    index_search = np.where(counts_im >= mincounts)
    if not len(index_search):
        print("NO PIXELS WITH ENOUGH COUNTS FOUND")
        return out_gain, out_counts, out_std

    for y, x in zip(*index_search):
        dist = data_cube[:, y, x]

        g, fit_g, success = fit_distribution(dist, pha_axis, mincounts = mincounts)

        if not success:
            continue

        # double-check
        if g.mean.value <= 3:
            sub_dist = dist - g(pha_axis) #g(np.arange(len(dist)))
            sub_dist[sub_dist < 0] = 0

            g2, fit2_g, success = fit_distribution(sub_dist, pha_axis, start_mean=15, mincounts = mincounts)

            if success and abs(g2.mean.value - g.mean.value) > 1:
                continue

        out_gain[y, x] = g.mean.value
        #out_counts[y, x] = dist.sum()
        out_std[y, x] = g.stddev.value

    if plot==True:
        plt.clf()
        plt.close()
        plt.plot(pha_axis, dist, 'k')
        plt.plot(pha_axis, g(pha_axis), 'r')
        plt.show()
        plt.clf()
        plt.close()

    return(out_gain, out_std)


def fit_distribution(dist, pha_axis,**kwargs):

    start_mean = kwargs.get('start_mean', None)
    start_amp = kwargs.get('start_amp', None)
    start_std = kwargs.get('start_std', None)
    mincounts = kwargs.get('mincounts', 10)
    
    x_vals =pha_axis #np.arange(len(dist))

    start_mean = start_mean or dist.argmax()
    start_amp = start_amp or int(max(dist))
    start_std = start_std or 1.05

    g_init = models.Gaussian1D(amplitude=start_amp,
                               mean=start_mean,
                               stddev=start_std,
                               bounds={'mean': [1, 30]})
    g_init.stddev.fixed = True

    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x_vals, dist)

    success = fit_ok(g, fit_g, start_mean, start_amp, start_std, mincounts = mincounts)

    return(g, fit_g, success)


def fit_ok(fit, fitter, start_mean, start_amp, start_std, mincounts = 10):

    # Check for success in the LevMarLSQ fitting
    if not fitter.fit_info['ierr'] in [1, 2, 3, 4]:
        return(False)

    # If the peak is too low
    if fit.amplitude.value < mincounts/np.sqrt(2.*math.pi)/fit.stddev.value:
        return(False)

    if not fit.stddev.value:
        return(False)

    # Check if fitting stayed at initial
    if not (start_mean - fit.mean.value):
        return(False)
    if not (start_amp - fit.amplitude.value):
        return(False)

    # Not sure this is possible, but checking anyway
    if np.isnan(fit.mean.value):
        return(False)
    if (fit.mean.value <= 0) or (fit.mean.value >= 31):
        return(False)

    return(True)


def create_per_mode_ccis_old(life_adj, hv, cenwave, fppos , segment , delta_t = 7,  compute_gain=False, start='2017-10-01'):
    
    if segment =='FUVB':
        ss = 'b'
    else:
        ss = 'a'

    ny = 1024
    nx = 16384

    #outfile
    outfile = '{}_lp{}_hv{}_c{}_fp{}_{}_dt{}'.format(start,life_adj, hv, cenwave, fppos, segment, delta_t)
        
    #Create the time axis
    cos_start = Time(start).mjd
    now = Time(datetime.datetime.now()).mjd
    time = np.arange(cos_start, now, delta_t) + delta_t/2
    ntime = len(time)

    #Create the arrays that will hold the counts, gain and exposure time

    pha_low = 1
    pha_high = 31
    npha = pha_high-pha_low + 2
    pha_axis = np.arange(pha_low, pha_high + 1, 1)

    counts = np.zeros((ntime, npha, ny, nx), dtype = 'float32')
    if compute_gain==True:
        gain = np.zeros((ntime, ny, nx),dtype = 'float32')
    exp_time = np.zeros(ntime, dtype = 'float32')

    #now get the observations with specified mode for each time interval
    
    for it in range(ntime):
        this_time = time[it]
        mint = this_time-delta_t/2
        maxt = this_time + delta_t/2

        if fppos ==1234:
            if segment =='FUVB':
                fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.life_adj==life_adj) & (FUVB_raw_headers.hvlevelb==hv) & (Observations.expstart >= mint) & (Observations.expstart < maxt))
            else:
                fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.life_adj==life_adj) & (FUVA_raw_headers.hvlevela==hv)& (Observations.expstart >= mint) & (Observations.expstart < maxt))
       
        else:
            if segment =='FUVB':
                fuv_obs = Observations.select().join(FUVB_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.fppos ==fppos) & (Observations.life_adj==life_adj) & (FUVB_raw_headers.hvlevelb==hv) & (Observations.expstart >= mint) & (Observations.expstart < maxt))
            else:
                fuv_obs = Observations.select().join(FUVA_raw_headers).where((Observations.detector=='FUV') & (Observations.cenwave==cenwave) & (Observations.fppos ==fppos) & (Observations.life_adj==life_adj) & (FUVA_raw_headers.hvlevela==hv)& (Observations.expstart >= mint) & (Observations.expstart < maxt))
    
        

        print("FOUND {} {}  Observations taken between {} and {} with cenwave {} FPPOS {} at lifetime postion {} and HV-{} = {}".format(len(fuv_obs), segment, mint, maxt, cenwave, fppos, life_adj, segment, hv))

        if len(fuv_obs)>0:

            for obs in fuv_obs:
                rootname = obs.rootname
                pid = int(obs.proposid)
                path = '/grp/hst/cos2/cosmo/{}/{}_corrtag_{}.fits.gz'.format(np.round(pid), rootname, ss)
                #For testing purposes, I am going to ignore proprietary data
                if os.access(path, os.R_OK):
                    corrtag = fits.open(path)
                    ctag = corrtag[1].data
                    cube = make_cube_pha(ctag, phamin = pha_low, phamax = pha_high)
                    counts[it,:,:,:] = counts[it,:,:,:] + cube
                    exp_time[it] = exp_time[it] + obs.exptime
                    
            if compute_gain==True:
                gain[it,:,:], gstd = measure_gainimage(cube[1:,:,:], phlow = pha_low, phhigh = pha_high)


   

    if (cenwave==1280) or (cenwave==800) or (cenwave==1105) or (cenwave==1105):
        opt_elem = 'G140L'
    if (cenwave ==1055) or (cenwave==1096) or (cenwave==1222) or (cenwave==1291) or (cenwave==1300) or (cenwave==1309) or (cenwave ==1318) or (cenwave==1327):
        opt_elem = 'G130M'
    if (cenwave ==1533) or (cenwave==1577) or (cenwave==1589) or (cenwave==1600) or (cenwave==1611) or (cenwave==1623):
        opt_elem = 'G160M'
            
    hduout = fits.HDUList(fits.PrimaryHDU())
    hduout[0].header['TELESCOP'] = 'HST'
    hduout[0].header['DETECTOR'] = 'FUV'
    hduout[0].header['OPT_ELEM'] = opt_elem
    hduout[0].header['CENWAVE'] = cenwave
    hduout[0].header['FPPOS']=fppos
    hduout[0].header['SEGMENT'] = segment
    hduout[0].header['HVLEVEL'] = hv
    hduout[0].header['FILETYPE'] = 'JRDCCI'
    hduout[0].header['EXPSTART'] = cos_start
    hduout[0].header['EXP_END'] = now
    hduout[0].header['DELTAT'] = delta_t
    hduout.append(fits.ImageHDU(time))
    hduout[1].header['EXTNAME']='TIME'
    hduout.append(fits.ImageHDU(pha_axis))
    hduout[2].header['EXTNAME'] = 'PHA'
    hduout.append( fits.ImageHDU(exp_time))
    hduout[3].header['EXTNAME'] = 'EXPTIME'
    hduout.append(fits.ImageHDU(counts))
    hduout[4].header['EXTNAME'] = 'COUNTS'
    hduout.writeto(outfile + '_CCI.fits', overwrite=True)
        
    if compute_gain==True:
        hdug = fits.HDUList(fits.PrimaryHDU())
        hdug[0].header['TELESCOP'] = 'HST'
        hdug[0].header['DETECTOR'] = 'FUV'
        hdug[0].header['OPT_ELEM'] = opt_elem
        hdug[0].header['CENWAVE'] = cenwave
        hdug[0].header['FPPOS']=fppos
        hdug[0].header['SEGMENT'] = segment
        hdug[0].header['HVLEVEL'] = hv
        hdug[0].header['FILETYPE'] = 'JRDCCI'
        hdug[0].header['EXPSTART'] = cos_start
        hdug[0].header['EXP_END'] = now
        hdug[0].header['DELTAT'] = delta_t 
        hdug.append(fits.ImageHDU(time))
        hdug[1].header['EXTNAME']='TIME'
        hdug.append(fits.ImageHDU(gain))
        hdug[2].header['EXTNAME'] = 'GAIN'
        hdug.writeto(outfile + '_GAIN.fits', overwrite=True)
            
                            
    if compute_gain==True:
        return(time, exp_time, counts, gain)
    else:
        return(time, exp_time, counts)
