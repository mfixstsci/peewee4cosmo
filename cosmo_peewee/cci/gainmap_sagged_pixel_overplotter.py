from astropy.io import fits
import os
import matplotlib.pyplot as plt
import glob

from ..database.models import get_settings

#-------------------------------------------------------------------------------
def make_overplot(gainsag_table, blue_modes=False):
    """
    Make plot of total gainmaps with gsagtab over plotted.

    Parameters
    ----------
    gainsag_table: str
        Path and filename of gainsag table
    blue_modes: bool
        Logic to make plots of blue modes or not.
    """
    
    #-- Get settings
    settings = get_settings()

    #-- Open GSAGTAB and gainmap.
    gsagtab = fits.open(gainsag_table)
    
    if not blue_modes:
        #-- Sort all of the FUV gainmaps by segment
        fuva_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 'CCI', '*_00_*gainmap*'))
        fuva_gainmaps.sort()
        fuvb_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 'CCI', '*_01_*gainmap*'))
        fuvb_gainmaps.sort()
        
        filename = 'gainmap_gsag_overplot.png'

    else:
        fuva_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 'CCI', '*_00_173*gainmap*'))
        fuva_gainmaps.sort()
        fuvb_gainmaps = glob.glob(os.path.join(settings['monitor_location'], 'CCI', '*_01_*gainmap*'))
        fuvb_gainmaps.sort()
        
        filename = 'gainmap_gsag_bluemodes_overplot.png'
     
    #-- Open the last created gainmaps for each segment.
    fuva_gainmap = fits.open(fuva_gainmaps[-1:][0])
    fuvb_gainmap = fits.open(fuvb_gainmaps[-1:][0])

    #-- Set the HVLVL and SEGMENT from the A/B gainmaps
    hv_lvl_a = int(fuva_gainmap[0].header['DETHV'])
    hv_lvl_b = int(fuvb_gainmap[0].header['DETHV'])
    
    #-- Use the HV from the last A/B gainmaps to open the total gainmaps.
    fuva_total_gainmap = fits.open(os.path.join(settings['monitor_location'], 'CCI','total_gain_{}.fits'.format(hv_lvl_a)))
    fuvb_total_gainmap = fits.open(os.path.join(settings['monitor_location'], 'CCI','total_gain_{}.fits'.format(hv_lvl_b)))

    #-- gsagtab keywords for high voltage are HVLEVEL[A/B].
    hvlvl_key = {'FUVA':'HVLEVELA',
                 'FUVB':'HVLEVELB'}

    #-- Find which extension in the GSAGTAB matches the HVLVL and SEGMENT for the gainmap.
    for ext in range(1, len(gsagtab)):
        gsagtab_segment = gsagtab[ext].header['segment']
        if (hv_lvl_a == gsagtab[ext].header[hvlvl_key[gsagtab_segment]]) and (gsagtab[ext].header['segment']=='FUVA'):
            fuva_gsag_ext = ext
        elif (hv_lvl_b == gsagtab[ext].header[hvlvl_key[gsagtab_segment]]) and (gsagtab[ext].header['segment']=='FUVB'):
            fuvb_gsag_ext = ext
        else:
            pass

    print(fuva_gsag_ext, fuvb_gsag_ext)
    #-- Create figure.
    f, axarr = plt.subplots(2, figsize=(20,15))

    #-- Set some plotting peramimeters.
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', lw=2)

    axes_font = 18
    title_font = 15

    #-- FUVA PANAL
    gm = axarr[0].imshow(fuva_total_gainmap['FUVALAST'].data, aspect='auto', cmap='binary')
    f.colorbar(gm, ax=axarr[0])
    axarr[0].scatter(gsagtab[fuva_gsag_ext].data['LX'], gsagtab[fuva_gsag_ext].data['LY'], marker='+', color='red', label='DQ 8192')
    axarr[0].legend(fontsize=15)

    axarr[0].set_title('SEGMENT: FUVA, HVLEVEL: {}, Gainmap: {}'.format(hv_lvl_a, os.path.basename(fuva_gainmaps[-1:][0])), fontsize=title_font, fontweight='bold')
    axarr[0].set_xlabel('X (Pixels)', fontsize=axes_font, fontweight='bold')
    axarr[0].set_ylabel('Y (Pixels)', fontsize=axes_font, fontweight='bold')

    axarr[0].set_xlim([400, 15500])
    axarr[0].set_ylim([200 , 800])
    #-- END FUVA PANEL

    #-- FUVB PANEL
    gm = axarr[1].imshow(fuva_total_gainmap['FUVBLAST'].data, aspect='auto', cmap='binary')
    f.colorbar(gm, ax=axarr[1])
    axarr[1].scatter(gsagtab[fuvb_gsag_ext].data['LX'], gsagtab[fuvb_gsag_ext].data['LY'], marker='+', color='red', label='DQ 8192')

    axarr[1].legend(fontsize=15)

    axarr[1].set_title('SEGMENT: FUVB, HVLEVEL: {}, Gainmap: {}'.format(hv_lvl_b, os.path.basename(fuvb_gainmaps[-1:][0])), fontsize=title_font, fontweight='bold')
    axarr[1].set_xlabel('X (Pixels)', fontsize=axes_font, fontweight='bold')
    axarr[1].set_ylabel('Y (Pixels)', fontsize=axes_font, fontweight='bold')

    axarr[1].set_xlim([400, 15400])
    axarr[1].set_ylim([350, 800])
    #-- END FUVB PANEL

    #-- Save figure.
    plt.savefig(os.path.join(settings['monitor_location'], 'CCI', filename))
    #-- Close.
    plt.close()
#-------------------------------------------------------------------------------