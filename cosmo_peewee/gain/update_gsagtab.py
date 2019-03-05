#!/usr/bin/env python

from __future__ import print_function, division

"""This module is meant to update a previously delivered gsagtab with new sagged pixels from gainmaps made on a two
 week cadence."""

__author__ = "Camellia Magness"
__email__ = "cmagness.stsci.edu"

import sys
import glob
import logging
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits

# global variables
OLD_GSAGTAB = "/user/cmagness/monitors/gainsag/23e16470l_gsag.fits"  # reference file date is 2018-03-13
GAINMAP_DIR = "/user/cmagness/monitors/gainsag/gainmaps/"
OUTDIR = "/user/cmagness/monitors/gainsag/"

logger = logging.getLogger(__gsagtab__)


# ----------------------------------------------------------------------------------------------------------------------


def files():
    """This function should be responsible for setting global variables or collecting necessary file system in some
    way."""
    # flags to consider
    # HV level, segment, blue mode/not
    parser = argparse.ArgumentParser(description="Script to update previously delivered gsagtab with new holes.")

    parser.add_argument("hv_level", type=str, help="hv level of gain to update")
    parser.add_argument("segment", type=str, help="detector segment of gain to update")
    parser.add_argument("--blue", action="store_true", default=False)

    args = parser.parse_args()

    all_gainmaps = glob.glob(GAINMAP_DIR + "*")
    if not args.blue:
        logger.info("COLLECTING GAINMAPS FOR HV LEVEL {} AND SEGMENT {}".format(args.hv_level, args.segment))
        gainmaps = [gainmap for gainmap in all_gainmaps if args.hv_level in gainmap if args.segment in gainmap]
    else:
        logger.info("BLUE MODES NOT SUPPORTED YET")
        logger.warning("GAINMAPS LIST IS EMPTY")
        gainmaps = []

    return gainmaps


# ----------------------------------------------------------------------------------------------------------------------


def stack_hdus(gainmaps):
    """This function stacks all binary maps from all relevant gainmaps into a dictionary with MJD values as keys."""
    logger.info("STACKING BINARY MAP HDUS")

    hdu_stack = {}
    for gainmap in gainmaps:
        with fits.open(gainmap) as f:
            prhd = f["PRIMARY"].header
            mjd = prhd["EXP_END"]
            binary_map = f["MASK DIF"].data
        hdu_stack[mjd] = binary_map

    logger.info("BINARY MAPS STACKED")

    return hdu_stack


# ----------------------------------------------------------------------------------------------------------------------


def find_sagged_pixels(hdu_stack, ymin, ymax, xmin, xmax):
    """This function should find the sagged pixels from the most recent gainmap that are not in the gsagtab."""
    logger.info("FINDING SAGGED PIXELS IN X REGION {}, {} AND Y REGION {}, {}".format(xmin, xmax, ymin, ymax))
    # most recent gainmap, by mjd
    mr_mjd = sorted(hdu_stack)[-1]
    search_region = hdu_stack[mr_mjd][ymin:ymax, xmin:xmax]
    sagged_pixels = np.nonzero(search_region)
    xpix = sagged_pixels[1]
    ypix = sagged_pixels[0]
    df_sagged = pd.DataFrame({"LX": xpix + xmin, "LY": ypix + ymin, "DATE": np.nan}, columns=["LX", "LY", "DATE"])
    # need the absolute value of the x, y pixels
    logger.info("FOUND {} SAGGED PIXELS NOT IN GAINSAG TABLE".format(len(df_sagged)))

    return df_sagged


# ----------------------------------------------------------------------------------------------------------------------


def find_sag_dates(df_sagged, hdu_stack, ymin, ymax, xmin, xmax):
    """This function should use the sagged pixels and scan the gainmaps to determine the date the pixel first sagged."""
    logger.info("FINDING DATES OF SAGGED PIXELS IN X REGION {}, {} AND Y REGION {}, {}".format(xmin, xmax, ymin, ymax))

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
        logger.info("FOUND PIXELS SAGGED AT DATE OF {}".format(mjd))
    logger.info("FOUND DATES FOR ALL SAGGED PIXELS")
    logger.info("ADDING DX, DY, DQ COLUMNS TO SAGGED PIXEL TABLE")

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
    logger.info("PLOTTING FOUND SAGGED PIXELS")
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
    logger.info("SAVING PLOT TO {}".format(OUTDIR + "gsag_regions.png"))
    plt.savefig(OUTDIR + "gsag_regions.png")


# ----------------------------------------------------------------------------------------------------------------------


def update_gsag_data(df_sagged):
    """This function should update the gsagtab with the newly sagged pixels and the date they were first sagged and
    create a new gsagtab in the OUTDIR."""

    # add conditional about whether segment is A or B -- HVLEVEL header keyword should be changed accordingly
    logger.info("FINDING CORRECT EXTENSION OF THE OLD GSAGTAB TO UPDATE")
    with fits.open(OLD_GSAGTAB) as f:

        # i am aware this is not... the best. suggestions?
        header_found = False
        extension = 1
        while not header_found:
            header = f[extension].header
            try:
                if (header["SEGMENT"] == "FUVB") & (header["HVLEVELB"] == 163):
                    logger.info("FOUND CORRECT EXTENSION")
                    header_found = True
            except:
                pass
            extension += 1
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
        f.writeto(OUTDIR + "new_gsagtab.fits", overwrite=True)
        logger.info("WROTE NEW GSAGTAB TO {}".format(OUTDIR + "new_gsagtab.fits"))


# ----------------------------------------------------------------------------------------------------------------------


def update_gsag_header():
    pass

# ----------------------------------------------------------------------------------------------------------------------


def main():
    # collect all gainmaps at hv level and segment of interest
    gainmaps = files()
    # stack all binary map extensions of these gainmaps
    hdu_stack = stack_hdus(gainmaps)
    # find pixels where holes exist in gainmap but not gsagtab & store coordinates
    df_sagged = find_sagged_pixels(hdu_stack, 475, 503, 1029, 14957)
    # ^^ these numbers are the lp4 region, need to add logic for defining different lps
    # compare this pixel list with each gainmap to find where it first appears, store date it first appears with pixel
    df_sagged_with_date = find_sag_dates(df_sagged, hdu_stack, 475, 503, 1029, 14957)
    plot_gainsag_regions(df_sagged_with_date)
    # add these pixels and date it first appears sagged to gsagtab. set dx and dy to 1
    update_gsag_data(df_sagged_with_date)

    # update history, do crds checks?
    update_gsag_header()

    # note for later: confirm that pixel is still sagged in the future and that it didn't show as sagged due to noise

    return 0


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    status = main()
    sys.exit(status)
