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
from astropy.io import fits


# global variables
OLD_GSAGTAB = "/user/cmagness/monitors/gainsag/23e1646sl_gsag.fits"  # reference file date is 2018-03-13
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


def find_sagged_pixels(hdu_stack, ymin, ymax):
    """This function should find the sagged pixels from the most recent gainmap that are not in the gsagtab."""
    logger.info("FINDING SAGGED PIXELS IN Y REGION {}, {}".format(ymin, ymax))
    # most recent gainmap, by mjd
    mr_mjd = sorted(hdu_stack)[-1]
    search_region = hdu_stack[mr_mjd][ymin:ymax, :]
    sagged_pixels = np.nonzero(search_region)
    xpix = sagged_pixels[1]
    ypix = sagged_pixels[0]
    df_sagged = pd.DataFrame({"xpix": xpix, "ypix": ypix + ymin}) # need the absolute value of the y pixels
    logger.info("FOUND {} SAGGED PIXELS NOT IN GAINSAG TABLE".format(len(df_sagged)))

    return df_sagged


# ----------------------------------------------------------------------------------------------------------------------


def find_sag_dates(df_sagged, hdu_stack, ymin, ymax):
    """This function should use the sagged pixels and scan the gainmaps to determine the date the pixel first sagged."""
    logger.info("FINDING DATES OF SAGGED PIXELS IN Y REGION {}, {}".format(ymin, ymax))

    for mjd, gain_hdu in sorted(hdu_stack.items()):
        search_region = gain_hdu[ymin: ymax, :]
        sagged_pixels = np.nonzero(search_region)
        xpix = sagged_pixels[1]
        ypix = sagged_pixels[0]
        df_comp = pd.DataFrame({"xpix": xpix, "ypix": ypix + ymin, "mjd_comp": mjd})
        df_sagged = pd.merge(df_sagged, df_comp, how="left", left_on=["xpix", "ypix"], right_on=["xpix", "ypix"])
        mask = (df_sagged["mjd"].isna()) & (df_sagged["mjd_comp"].notnull())
        df_sagged["mjd"].mask(mask, df_sagged["mjd_comp"], inplace=True)
        df_sagged = df_sagged.drop(columns=["mjd_comp"])
        logger.info("FOUND PIXELS SAGGED AT DATE OF {}".format(mjd))
    logger.info("FOUND DATES FOR ALL SAGGED PIXELS")

    return df_sagged


# ----------------------------------------------------------------------------------------------------------------------


def update_gsag(full_pixel_table):
    """This function should update the gsagtab with the newly sagged pixels and the date they were first sagged and
    create a new gsagtab in the OUTDIR."""
    pass


# ----------------------------------------------------------------------------------------------------------------------


def main():
    pass
    # collect all gainmaps at hv level and segment of interest
    gainmaps = files()
    # stack all binary map extensions of these gainmaps
    hdu_stack = stack_hdus(gainmaps)
    # find pixels where holes exist in gainmap but not gsagtab & store coordinates
    df_sagged = find_sagged_pixels(hdu_stack, 475, 503)
    # ^^ these numbers are the lp4 region, need to add logic for defining different lps
    # compare this pixel list with each gainmap to find where it first appears, store date it first appears with pixel
    df_sagged_with_date = find_sag_dates(df_sagged, hdu_stack, 475, 503)

    # note for later: confirm that pixel is still sagged in the future and that it didn't show as sagged due to noise

    # add these pixels and date it first appears sagged to gsagtab. set dx and dy to 1
    update_gsag(df_sagged_with_date)
    return 0


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    status = main()
    sys.exit(status)

