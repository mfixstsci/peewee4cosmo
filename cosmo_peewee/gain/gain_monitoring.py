from astropy.io import fits, ascii
from matplotlib import pyplot as plt
from astropy.convolution import convolve, Box1DKernel
import numpy as np
#import scipy.misc
#import scipy.ndimage
import skimage.transform
import glob
from astropy.time import Time
import matplotlib
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 
import sys
sys.path.append('/user/duval/COS/common')
sys.path.append("PEEWEE4COSMO/cosmo_peewee/database/")
from cos_utils import *
from models import Gain, Gain_Trends, Observations, Flagged_Pixels, get_database
import datetime

def check_gain_table(x,y):
    db = get_database()
    db.connect()
    glist = Gain.select().distinct().where((Gain.y==y) & (Gain.x==x) & (Gain.segment=='FUVB') & (Gain.hv_lvl ==163))
    #glist = Gain.select().distinct().where((Gain.gain>=1) & (Gain.gain <=3))

    #glist = Gain.select().limit(100)
    #print(glist)

    nobs = len(glist)
    expstarts = np.zeros(nobs,dtype = 'float32')
    gains = np.zeros(nobs, dtype = 'float32')
    counts = np.zeros(nobs, dtype = 'float32')

    pos = 0
    
    for g in glist:
        print("")
        print("X", g.x)
        print("Y ", g.y)
        print("Gain ", g.gain)
        print("Counts ", g.counts)
        print("HV ", g.hv_lvl)
        print("SEG ", g.segment)
        print("EXPSTART ", g.expstart)
        expstarts[pos] = g.expstart
        gains[pos] = g.gain
        counts[pos] = g.counts
        pos = pos + 1


    im = fits.open('/grp/hst/cos/Monitors/CCI/proj_bad_FUVB_163.fits')
    slopes = im['SLOPE'].data
    intercepts = im['INTERCEPT'].data
    projbad = im['PROJBAD'].data

    plt.clf()
    plt.close()
    plt.plot(expstarts, gains, 'ko')
    time = np.arange(min(expstarts), max(expstarts), 1)
    plt.plot(time, slopes[y*2, x*8]*time + intercepts[y*2, x*8], 'r-')
    plt.show()
    plt.clf()
    plt.close()
        

    db.close()

def check_gain_trend_table(x,y):
    db = get_database()
    db.connect()
    #glist = Gain.select().distinct().where((Gain.y==y) & (Gain.x==x) & (Gain.segment=='FUVB') & (Gain.hv_lvl ==163))
    #glist = Gain_Trends.select().distinct().where((Gain_Trends.segment=='FUVB') & (Gain_Trends.hv_lvl == 163) & (Gain_Trends.x==x) & (Gain_Trends.y == y))
    glist = Gain_Trends.select().where((Gain_Trends.segment=='FUVB'))

    hvs = []

    #glist = Gain.select().limit(100)
    #print(glist)


    for g in glist:
        #print("")
        #print(g.x, g.y)
        #print(g.hv_lvl)
        #print(g.slope, g.intercept)
        #print(g.proj_bad_mjd)
        hvs.append(g.hv_lvl)

    db.close()

    return(hvs)

    
    
    

def get_cos_hv_history(segment, this_time = None):

    now  = datetime.datetime.now()
    now_string = str(now.year) + '-'  + str(now.month) + '-' + str(now.day)

    if segment =='FUVA':

        tmin = ['2009-08-12', '2012-03-26', '2012-07-23' , '2014-11-03', '2015-02-09', '2017-10-02']
        tmax = ['2012-03-26', '2012-07-23', '2014-11-03', '2015-02-09', '2017-10-02', now_string]
        hv  = [169, 178, 167, 173, 167, 163]

    else:

        tmin = ['2009-08-12', '2011-08-03', '2012-07-23', '2013-06-24', '2014-07-21', '2015-02-09', '2016-01-18', '2016-10-17', '2017-10-02']
        tmax = ['2011-08-03', '2012-07-23', '2013-06-24', '2014-07-21', '2015-02-09', '2016-01-18', '2016-10-17', '2017-10-02', now_string]
        hv= [167, 175, 163, 169, 175, 163, 169 , 175, 163]

    if this_time !=None:

        mintimes = Time(tmin, format= 'iso').mjd
        maxtimes = Time(tmax, format = 'iso').mjd
        get_this_time = Time(this_time, format = 'iso')
        get_this_time = get_this_time.mjd

        dif = get_this_time-mintimes
        dif[np.where(dif <0)] = 1.e10
        index = np.argmin(dif)
        this_hv = hv[index]

        return(tmin, tmax, hv, this_hv)
    else:
        return(tmin, tmax, hv)
        

    

def read_gsag_table(segment, hv, mindate='2009-01-01', file = '', region_file = '', color='blue'):

    """
    This program reads the gain sag table adn outputs the list of gain-sag holes (lx,ly, dx, dy) that appears on or after mindate
    
    INPUTS:
    
    segment ("FUVA" or "FUVB")
    hv (e.g., 163)
    
    file is the GSAG reference file
    region_file is hte name of an output region file which as the boxes coresponding to gain sag holes that have appeared after mindate.
    """

    start = Time(mindate).mjd
    
    if segment=='FUVB':
        ss = 'B'
    else:
        ss = 'A'

    if file != '':
        t = fits.open(file)
        print("READING ", file)
    else:
        t = fits.open("/grp/crds/hst/references/hst/23e16470l_gsag.fits")
        print("READING /grp/crds/hst/references/hst/23e16470l_gsag.fits")

    #find the right extension

    for i in range(1,len(t),1):
        
        seg = t[i].header['SEGMENT']
        if seg==segment:
            this_hv = t[i].header['HVLEVEL' + ss]
            if hv==this_hv:
                extension =i
                break
    gsag_tab = t[extension].data

    dates = Time(gsag_tab['DATE'], format = 'mjd').mjd

    index =np.where(dates >= start)
    index = index[0]
    lx = gsag_tab['LX'][index]
    ly = gsag_tab['LY'][index]
    dx = gsag_tab['DX'][index]
    dy = gsag_tab['DY'][index]
    date = dates[index]

    ngsag =  len(lx)
    
    if region_file != '':
        f = open(region_file, 'w')
        for i in range(ngsag):
            cx = lx[i] + dx[i]/2 
            cy = ly[i] + dy[i]/2
            reg_string = 'box({}, {}, {}, {}, 0) #color={}'.format(cx, cy, dx[i], dy[i], color)
            f.write(reg_string+'\n')
            
        f.close

    return(lx, ly, dx, dy, date)
    
    
            

def get_yrange(lp, segment):
    if lp == 1:
        if segment =='FUVB':
            ymin = 531
            ymax = 561
        if segment =='FUVA':
            ymin = 472
            ymax= 502

    if lp == 2:
        if segment =='FUVB':
            ymin = 572
            ymax= 595
        if segment =='FUVA':
            ymin = 507
            ymax = 550

    if lp == 3:
        if segment =='FUVB':
            ymin = 505
            ymax = 530
        if segment =='FUVA':
            ymin  = 447
            ymax = 472

    if lp==4:
        if segment=='FUVB':
            ymin = 470
            ymax = 493
        if segment=='FUVA':
            ymin = 406
            ymax = 441

    return(ymin, ymax)


def check_ccis_for_data(year, lp, segment, hv):

    if segment == 'FUVA':
        cci_type = '00'
    else:
        cci_type = '01'
        
    files = glob.glob('/grp/hst/cos/Monitors/CCI/l_' + '{}'.format(year) + '*_' + cci_type + '_' + '{}'.format(hv) + '_cci_gainmap.fits.gz')
    doys = []


    ymin, ymax  = get_yrange(lp, segment)
    ymin = ymin/2
    ymax = ymax/2 #to account for 512 sampling grid in CCIs
    for file in files:
        g = fits.open(file)
        c = g['COUNTS'].data
        lpc = np.nansum(c[ymin:ymax,:])

        if lpc > 5e-5*16384*(ymax-ymin)*g[0].header['EXPTIME']:
            splt = file.split('/')
            root = splt[len(splt)-1]
            this_doy = root[6:9]
            this_doy =  np.float(this_doy)
            doys.append(this_doy)
            print('COUNTS ', this_doy, lpc)

    return(doys)


def make_colorful_plot_lp3(year, lp, segment, hv):

    #doys = check_ccis_for_data(2018, 3, segment, hv)

    doys = [93, 274]
    print(doys)
    hvs = [hv,]*len(doys)
    print(hvs)
    
    years = [2018,]*len(doys)
    print(years)
    
    colorful_plot(3, segment, hvs, years, doys, outfile = 'colorful_plot_' + segment + '_lp3_oct2018.pdf')
            
            


def get_lp_boundaries_1dx(lp, segment):

    dir = '/grp/crds/hst/references/hst/'
    lp4_1dx = '23n1744ol_1dx.fits'
    lp3_1dx = 'z2d19237l_1dx.fits'

    if lp == 3:
        dx_file = lp3_1dx
    else:
        dx_file = lp3_1dx

    t = fits.open(dir + dx_file)
    t = t[1].data
    index = np.where((t['APERTURE']=='PSA') & (t['SEGMENT']==segment))
    index =index[0]
    ymin = 2000
    ymax = 0
    for i in range(len(index)):
        ylow = t['B_SPEC'][index[i]] + 16384.*t['SLOPE'][index[i]] - t['HEIGHT'][index[i]]/2
        yhigh = t['B_SPEC'][index[i]] + 16384.*t['SLOPE'][index[i]] + t['HEIGHT'][index[i]]/2
        if np.min(ylow)<ymin:
            ymin = np.min(ylow)
        if np.max(yhigh)>ymax:
            ymax  = np.max(yhigh)

    return(ymin, ymax)


def plot_gain_profiles(segment, hv, cenwave, shifty = 0):


    gain_file = "/grp/hst/cos/Monitors/CCI/total_gain_" + hv + ".fits"
    gain = fits.open(gain_file)
    gain = gain[segment + 'LAST'].data


    if segment =='FUVB':
        ymin = 460
        ymax = 530
    else:
        ymin = 400
        ymax = 470

    if shifty ==0:
        key = ''
    else:
        key = '_shift{}'.format(shifty)
        
    plot_spectral_image(gain, 'plot_gain_' + segment + '_' + hv + '_prof_'+ str(cenwave) + key + '.pdf', profile  =True, cenwave = cenwave, segment = segment, ymin = ymin, ymax = ymax, shifty = shifty)


def map_gain(segment, hv):

    gain_file = 'total_gain_' + '{}'.format(hv) + '.fits'
    gain=fits.open(gain_file)
    gain = gain[segment + 'LAST'].data

    ymin3, ymax3 = get_lp_boundary_1dx(3, segment)
    ymin4, ymax4 = get_lp_boundary_1dx(4, segment)
    
    plt.clf()
    plt.close()
    fig = plt.figure(figsize=(15, 8))
    plt.imshow(gain, vmin = 0, vmax = 12, cmap = 'viridis_r', aspect = 'auto', origin = 'lower')
    plt.axhline(ymin4, color = 'red')
    plt.axhline(ymax4, color = 'red')
    plt.axhline(ymin3, color = 'magenta')
    plt.axhline(ymax3, color = 'magenta')
    
    plt.xlabel("XCORR", fontsize = 18)
    plt.ylabel("YCORR", fontsize = 18)
    cbar = plt.colorbar()
    cbar.set_title("MODAL GAIN", fontsize = 18)
    plt.savefig("map_gain_" + segment + "_" + '{}'.format(hv) + ".pdf", format = 'pdf', dpi = 1000)
    plt.clf()
    plt.close()


def resample_gain(gain_file,extension, data =None):

    if data != None:
        gain = data
    else:
        g = fits.open(gain_file)
        gain = g[extension].data
    gain = skimage.transform.resize(gain, (1024, 16384), preserve_range=True, mode = 'constant', cval = 0., order = 0)
    return(gain)

def get_gain_file(segment, hv, year, doy, from_total=True):

    dir = '/grp/hst/cos/Monitors/CCI/'
    if year == 'latest':
        gain_files = 'total_gain_' + '{}'.format(hv) + '.fits'
       
    else:

        ndates = len(doy)
        if segment == 'FUVA':
            cci_type = '00'
        else:
            cci_type = '01'

        gain_files = []
        actual_doy = []
        mjd = []
        
        for i in range(ndates):

            if from_total==False:

                files = glob.glob(dir + "l_" + '{}'.format(year[i]) + '*_' + cci_type + '_' + '{}'.format(hv[i]) + '_cci_gainmap.fits.gz')
                doys = np.zeros(len(files),dtype = 'uint32')
                for j in range(len(files)):
                    file = files[j]
                    splt = file.split('/')
                    root = splt[len(splt)-1]
                    this_doy = root[6:9]
                    doys[j] = this_doy

            else:
            
                files = glob.glob("/grp/hst/cos/duval/TOTAL_GAIN/SCALED_TOTAL_GAIN_FOR_MOVIE/total_gain_{}_{}_{}-*-*_scaled.fits".format(segment, hv[i], year[i]))
                doys  = np.zeros(len(files), dtype= 'float32')
                for j in range(len(files)):
                    file = files[j]
                    splt = file.split('/')
                    root = splt[len(splt)-1]
                    date = root.split('_')
                    date = date[len(date)-2]
                    date_time = Time(date, format = 'iso')
                    date_time.format = 'decimalyear'
                    this_doy= (date_time.value-np.floor(date_time.value))*365. + 1.
                    doys[j] = this_doy

            
            dif = np.abs(doy[i]-doys)
            this_doy = doys[np.argmin(dif)]
            actual_doy.append(this_doy)
            gain_file = files[np.argmin(dif)]
            g = fits.open(gain_file)
            if from_total==False:
                this_mjd = g[0].header['EXPSTART']
            else:
                this_mjd = g[0].header['EXP_END']
            mjd.append(this_mjd)
            gain_files.append(gain_file)
            
    return(gain_files, actual_doy, mjd)



def make_plot_for_dave():
    lp= 4
    segment = 'FUVB'
    hvs = [163, 163, 163]
    years = [2018, 2018, 2018]
    doys = [210, 240, 260]
    colorful_plot(lp, segment, hvs, years, doys, outfile = 'dave_plot.pdf')
            
def colorful_plot(lp, segment, hvs, years, doys,outfile = '', smooth = 51, from_total=True):

    #Need for loop over all dates
    n = len(years)

    plt.clf()
    plt.close()

    #fig,ax = plt.subplots(2, sharex = True,figsize = (10, 16))
    fig = plt.figure(figsize = (10,8))
    gain_files, actual_doys, mjds = get_gain_file(segment, hvs, years, doys, from_total = from_total)
    ymin, ymax = get_yrange(lp, segment)

    print("GAIN FILES ", gain_files)
    print("DOY ", actual_doys)

    ycorr = np.arange(1024)
    xcorr = np.arange(16384)
    
    yy = ycorr[ymin:ymax]

    
    for iy in range(len(years)):
        gain_file = gain_files[iy]
        if from_total==False:
            gain = resample_gain(gain_file, 'MOD_GAIN')
            counts = resample_gain(gain_file, 'COUNTS')
        else:
            g = fits.open(gain_file)
            gain = g['TOTAL GAIN'].data
            counts = g['TOTAL COUNTS'].data
        #if iy == 0:
        #    hdu = fits.PrimaryHDU(gain)
        #    hdu.writeto("test_gain.fits", overwrite=True)
        
        yslice = np.zeros(16384, dtype = 'float32')
        gslice = np.zeros(16384, dtype = 'float32')
        cslice = np.zeros(16384, dtype = 'float32')

        at = Time(mjds[iy], format = 'mjd')
        date = '{}'.format(at.datetime.month) + '/' +'{}'.format(at.datetime.day) + '/' + '{}'.format(at.datetime.year)
        
        for ix in range(len(xcorr)):
            x = xcorr[ix]
            gcol = gain[ymin:ymax, x]
            ccol = counts[ymin:ymax, x]
            valid  = np.where(gcol>0)
            valid = valid[0]
            if len(valid)> 0:
                yslice[ix] =yy[np.argmin(gcol[valid])]
                gslice[ix] = gain[yslice[ix], x]
                cslice[ix] = counts[yslice[ix], x]

        good = np.where(gslice > 0.)
        good = good[0]
        plot_gslice = gslice[good]
        
        if smooth >1:
            plot_gslice  = convolve(plot_gslice, Box1DKernel(smooth))
        
            
        plt.plot(xcorr[good], plot_gslice, '-', label = date + ' HV'  + '{}'.format(hvs[iy]), alpha = 0.8)
        #ax[1].plot(xcorr[good], cslice[good], '-', label = date + ' HV'  + '{}'.format(hvs[iy]))

    #ax[0].set_xticklabels([])
    #ax[1].set_xlabel("XCORR", fontsize = 20)
    #ax[0].set_ylabel("MODAL GAIN", fontsize = 20)
    #ax[1].set_ylabel("COUNTS", fontsize = 20)
    #fig.subplots_adjust(left = 0.1, right = 0.95, bottom = 0.1, top = 0.95, hspace = 0)
    #ax[0].legend(fontsize = 16, loc = 'upper left')

    #plt.xticklabels([])
    plt.xlabel("XCORR", fontsize = 20)
    plt.ylabel("MODAL GAIN", fontsize = 20)
    #ax[1].set_ylabel("COUNTS", fontsize = 20)
    #fig.subplots_adjust(left = 0.1, right = 0.95, bottom = 0.1, top = 0.95, hspace = 0)
    plt.legend(fontsize = 16, loc = 'upper left')
    plt.tight_layout()

    if outfile != '':
        plt.savefig(outfile, format = 'pdf', dpi = 1000)
    else:
        plt.show()
    
    plt.clf()
    plt.close()
