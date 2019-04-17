#!/usr/bin/env python

from __future__ import print_function, division

"""This module is meant to update a previously delivered gsagtab with new sagged pixels from gainmaps made on a two
 week cadence."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import os
import sys
import glob
import logging
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from datetime import datetime
from shutil import copy

# global variables
OLD_GSAGTAB = "/user/cmagness/monitors/gainsag/23e16470l_gsag.fits"  # reference file date is 2018-03-13
GAINMAP_DIR = "/user/cmagness/monitors/gainsag/gainmaps/"
OUTDIR = "/user/cmagness/monitors/gainsag/"


# ----------------------------------------------------------------------------------------------------------------------


class Info:
    """This class is for storing information to be passed around by this module from the initialization."""

    # add some support for optional flag (like blue mode?)

    def __init__(self, hv_level, segment, v_calcos, lp, merger):
        self.hv = hv_level
        self.segment = segment
        self.calcos = v_calcos
        self.merger = merger
        self.gainmaps = []

        self.xmin = 1029
        self.xmax = 14957

        if lp == 4:
            self.ymin = 475
            self.ymax = 503
        else:
            print("That lifetime position isn't technically supported yet. Setting to default LP4 values")
            self.ymin = 475
            self.ymax = 503

        # there has to be a better way to store this because this just get updated by hand right now
        self.pedigree = 'INFLIGHT 06/08/2009 24/02/2019'
        self.descrip = 'Updated to add new gainsagged regions at LP4 for HV Level of 163'
        self.history = ['--------------------- March 2019 --------------------------',
                        'These new regions were derived by C. Magness and J. Roman-Duval.']

    # tbh i made this a method largely because i was like a class that just initializes stuff is kind of boring
    def add_gainmaps(self, gainmaps):
        self.gainmaps = gainmaps


# ----------------------------------------------------------------------------------------------------------------------


def files():
    """This function should be responsible for setting global variables or collecting necessary file system in some
    way."""

    parser = argparse.ArgumentParser(description="Script to update previously delivered gsagtab with new holes.")

    parser.add_argument("hv_level", type=str, help="hv level of gain to update")
    parser.add_argument("segment", type=str, help="detector segment of gain to update")
    parser.add_argument("v_calcos", type=str, help="most recent version of calcos")
    parser.add_argument("lifetime_pos", type=int, help="LP to update, accepts 1-4")
    parser.add_argument("--blue", action="store_true", default=False, help="optional flag for blue mode gsagtab")
    # ^ this doesn't do anything at the moment but whatever

    args = parser.parse_args()

    merger = input("Please enter your name (first initial and last, ie. C. Magness): ")

    info = Info(args.hv_level, args.segment, args.v_calcos, args.lifetime_pos, merger)

    all_gainmaps = glob.glob(GAINMAP_DIR + "*")
    if not args.blue:
        logging.info("COLLECTING GAINMAPS FOR HV LEVEL {} AND SEGMENT {}".format(args.hv_level, args.segment))
        gainmaps = [gainmap for gainmap in all_gainmaps if args.hv_level in gainmap if args.segment in gainmap]
    else:
        logging.info("BLUE MODES NOT SUPPORTED YET")
        logging.warning("GAINMAPS LIST IS EMPTY")
        gainmaps = []

    info.add_gainmaps(gainmaps)

    return info


# ----------------------------------------------------------------------------------------------------------------------


def stack_hdus(gainmaps):
    """This function stacks all binary maps from all relevant gainmaps into a dictionary with MJD values as keys."""

    logging.info("STACKING BINARY MAP HDUS")

    hdu_stack = {}
    for gainmap in gainmaps:
        with fits.open(gainmap) as f:
            prhd = f["PRIMARY"].header
            mjd = prhd["EXP_END"]
            binary_map = f["MASK DIF"].data
        hdu_stack[mjd] = binary_map

    logging.info("BINARY MAPS STACKED")

    return hdu_stack


# ----------------------------------------------------------------------------------------------------------------------


def find_sagged_pixels(hdu_stack, ymin, ymax, xmin, xmax):
    """This function should find the sagged pixels from the most recent gainmap that are not in the gsagtab."""

    logging.info("FINDING SAGGED PIXELS IN X REGION {}, {} AND Y REGION {}, {}".format(xmin, xmax, ymin, ymax))

    # most recent gainmap, by mjd
    mr_mjd = sorted(hdu_stack)[-1]
    search_region = hdu_stack[mr_mjd][ymin:ymax, xmin:xmax]
    sagged_pixels = np.nonzero(search_region)
    xpix = sagged_pixels[1]
    ypix = sagged_pixels[0]
    df_sagged = pd.DataFrame({"LX": xpix + xmin, "LY": ypix + ymin, "DATE": np.nan}, columns=["LX", "LY", "DATE"])
    # need the absolute value of the x, y pixels
    logging.info("FOUND {} SAGGED PIXELS NOT IN GAINSAG TABLE".format(len(df_sagged)))

    return df_sagged


# ----------------------------------------------------------------------------------------------------------------------


def find_sag_dates(df_sagged, hdu_stack, ymin, ymax, xmin, xmax):
    """This function should use the sagged pixels and scan the gainmaps to determine the date the pixel first sagged."""

    logging.info("FINDING DATES OF SAGGED PIXELS IN X REGION {}, {} AND Y REGION {}, {}".format(xmin, xmax, ymin, ymax))

    for mjd, gain_hdu in sorted(hdu_stack.items()):
        search_region = gain_hdu[ymin: ymax, xmin:xmax]
        sagged_pixels = np.nonzero(search_region)
        xpix = sagged_pixels[1]
        ypix = sagged_pixels[0]
        df_comp = pd.DataFrame({"LX": xpix + xmin, "LY": ypix + ymin, "DATE_COMP": mjd}, columns=["LX", "LY",
                                                                                                  "DATE_COMP"])
        df_sagged = pd.merge(df_sagged, df_comp, how="left", left_on=["LX", "LY"], right_on=["LX", "LY"])
        mask = (df_sagged["DATE"].isna()) & (df_sagged["DATE_COMP"].notnull())
        df_sagged["DATE"].mask(mask, df_sagged["DATE_COMP"], inplace=True)
        df_sagged = df_sagged.drop(columns=["DATE_COMP"])
        logging.info("FOUND PIXELS SAGGED AT DATE OF {}".format(mjd))
    logging.info("FOUND DATES FOR ALL SAGGED PIXELS")
    logging.info("ADDING DX, DY, DQ COLUMNS TO SAGGED PIXEL TABLE")

    df_sagged["DX"] = 1
    df_sagged["DY"] = 1
    df_sagged["DQ"] = 8192
    columns = ["DATE", "LX", "LY", "DX", "DY", "DQ"]
    df_sagged = df_sagged[columns]
    df_sagged = df_sagged.apply(pd.to_numeric)

    return df_sagged


# ----------------------------------------------------------------------------------------------------------------------


def plot_gainsag_regions(df_sagged, ymin, ymax, xmin, xmax):
    """This function should plot the gainsagged regions found in the find_sagged_pixels and find_sag_dates functions."""

    logging.info("PLOTTING FOUND SAGGED PIXELS")

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    plt.title("Gainsag Regions in Gainmaps Missing in GSAGTAB")

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin - 5, ymax + 5)

    major_ticks_x = np.arange(xmin, xmax, 1000)
    major_ticks_y = np.arange(ymin, ymax, 5)
    minor_ticks_x = np.arange(xmin, xmax, 100)
    minor_ticks_y = np.arange(ymin, ymax, 1)

    ax.set_xticks(major_ticks_x)
    ax.set_xticks(minor_ticks_x, minor=True)
    ax.set_yticks(major_ticks_y)
    ax.set_yticks(minor_ticks_y, minor=True)

    ax.grid(which="both")
    ax.grid(which="major", linestyle="--", alpha=0.5)
    ax.grid(which="minor", linestyle="--", alpha=0.2)

    im = ax.scatter(df_sagged["LX"], df_sagged["LY"], s=1, c=df_sagged["DATE"], cmap="YlOrRd")
    cbar = fig.colorbar(im)
    cbar.set_label("DATE (MJD)")
    plt.show()
    logging.info("SAVING PLOT TO {}".format(OUTDIR + "gsag_regions.png"))
    plt.savefig(OUTDIR + "gsag_regions.png")


# ----------------------------------------------------------------------------------------------------------------------


def update_gsag_data(df_sagged, info):
    """This function should update the gsagtab with the newly sagged pixels and the date they were first sagged and
    create a new gsagtab in the OUTDIR. It should update the header by calling update_gsag_header."""

    # add conditional about whether segment is A or B -- HVLEVEL header keyword should be changed accordingly
    logging.info("FINDING CORRECT EXTENSION OF THE OLD GSAGTAB TO UPDATE")

    duplicate = copy(OLD_GSAGTAB, "./duplicate.fits")
    # ^ this has to be done so that running this multiple times you aren't overwriting the original file. there's
    # definitely a better way to do this but this was a quick fix. duplicate is deleted later after file is written

    with fits.open(duplicate) as f:

        # i am aware this is not... the best. suggestions?
        header_found = False
        extension = 0
        while not header_found:
            header = f[extension].header
            try:
                if (header["SEGMENT"] == "FUVB") & (header["HVLEVELB"] == 163):
                    logging.info("FOUND CORRECT EXTENSION")
                    header_found = True
                    break
            except:
                pass
            extension += 1
            # I KNOW THIS IS BAD
        data = f[extension].data

        date_col = fits.Column(name="DATE", format="D", unit="MJD", array=np.append(data["DATE"], df_sagged["DATE"]))
        lx_col = fits.Column(name="LX", format="J", unit="pixel", array=np.append(data["LX"], df_sagged["LX"]))
        ly_col = fits.Column(name="LY", format="J", unit="pixel", array=np.append(data["LY"], df_sagged["LY"]))
        dx_col = fits.Column(name="DX", format="J", unit="pixel", array=np.append(data["DX"], df_sagged["DX"]))
        dy_col = fits.Column(name="DY", format="J", unit="pixel", array=np.append(data["DY"], df_sagged["DY"]))
        dq_col = fits.Column(name="DQ", format="J", array=np.append(data["DQ"], df_sagged["DQ"]))
        cols = fits.ColDefs([date_col, lx_col, ly_col, dx_col, dy_col, dq_col])
        bintable = fits.BinTableHDU.from_columns(cols)

        f[extension].data = bintable.data
        bitpix = fits.getval(duplicate, 'BITPIX')
        # update header in place here before writing out the new file
        update_gsag_header(f["PRIMARY"], info, bitpix)
        f.writeto(OUTDIR + "new_gsagtab.fits", overwrite=True, checksum=True)
        os.remove(duplicate) # BAD

    logging.info("WROTE NEW GSAGTAB TO {}".format(OUTDIR + "new_gsagtab.fits"))


# ----------------------------------------------------------------------------------------------------------------------


def update_gsag_header(hdu, info, bitpix):
    """This function updates the header of an hdu in place. Should be the primary header. Adapted from another reference
    file editing module written by E. Frazer, C. Magness, and R. Plesha."""

    # appending a new description of the work that has been done to the existing header
    # also updating the pedigree, vcalcos, comment, and history
    # each history has to be less than 72 characters

    logging.info("UPDATING GSAGTAB HEADER")

    def descrip_keyword(keyword_text):
        required_length = 67
        current_length = len(keyword_text)
        remaining = required_length - current_length

        if remaining < 0:
            raise ValueError('DESCRIP keyword is too long (maximum is {}). Current length: {}'.format(required_length,
                                                                                                      current_length))
        elif remaining == 0:
            return keyword_text
        else:
            return keyword_text + '-' * remaining

    checked_descrip = descrip_keyword(info.descrip)

    del hdu.header['COMMENT']
    del hdu.header['FILENAME']

    for entry in info.history:
        hdu.header['HISTORY'] = entry

    try:
        del hdu.header['GIT_TAG']
    except KeyError:
        # Some reference files don't include a GIT_TAG yet, but if the keyword
        #   already exists it will place 2 in there which is unwanted.
        pass

    date_time = datetime.now().replace(microsecond=0).isoformat()
    hdu.header.set('DATE', date_time, 'Creation UTC (CCCC-MM-DD) date of FITS header')
    hdu.header.set('DESCRIP', value=checked_descrip)  # Reason for delivery ; needs to be 67 characters long exactly
    hdu.header.set('PEDIGREE', value=info.pedigree)
    hdu.header.set('BITPIX', value=bitpix)
    hdu.header['VCALCOS'] = info.calcos  # minimum version of calcos that can be run with this file
    hdu.header.add_comment('Reference file created by {}'.format(info.merger), before="HISTORY")

    logging.info("GSAGTAB HEADER UPDATED")

# ----------------------------------------------------------------------------------------------------------------------


def main():

    logging.basicConfig(filename='gsagtab.log', level=logging.DEBUG)
    # collect all gainmaps at hv level and segment of interest, make an object to hold that and relevant header info
    info = files()
    # stack all binary map extensions of these gainmaps
    hdu_stack = stack_hdus(info.gainmaps)
    # find pixels where holes exist in gainmap but not gsagtab & store coordinates
    df_sagged = find_sagged_pixels(hdu_stack, info.ymin, info.ymax, info.xmin, info.xmax)
    # ^^ these numbers are the lp4 region, need to add logic for defining different lps
    # compare this pixel list with each gainmap to find where it first appears, store date it first appears with pixel
    df_sagged_with_date = find_sag_dates(df_sagged, hdu_stack, info.ymin, info.ymax, info.xmin, info.xmax)
    plot_gainsag_regions(df_sagged_with_date, info.ymin, info.ymax, info.xmin, info.xmax)
    # add these pixels and date it first appears sagged to gsagtab. set dx and dy to 1
    # update header within this function and generate a new gsagtab
    update_gsag_data(df_sagged_with_date, info)

    # note for later: confirm that pixel is still sagged in the future and that it didn't show as sagged due to noise

    return 0


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    status = main()
    sys.exit(status)
