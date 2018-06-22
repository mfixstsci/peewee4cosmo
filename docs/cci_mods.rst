CCI Monitor
===========
The Cumulative Counts Image (CCI) monitor tracks CCI images. These images are roughly a weeks worth
of exposure for COS and are an image of counts in xcorr, ycorr space. This monitor 
creates modal gain maps for each CCI that is created and stores every modal gain measurement in a 
COSMO database table. The modal gain maps are a crucial part of tracking any pixels that drop below 
the modal gain threshold. These pixels are flagged and stored in a COSMO database table that is then used 
to create the gain sag reference table (GSAGTAB) used by the pipeline to properly flag sagging pixels.

find_bad_pix
------------
.. automodule:: cosmo_peewee.cci.find_bad_pix
   :members:

gainmap
-------
.. automodule:: cosmo_peewee.cci.gainmap
   :members:

gainsag
-------
.. automodule:: cosmo_peewee.cci.gainsag
   :members:

gainmap_sagged_pixel_overplotter
--------------------------------
.. automodule:: cosmo_peewee.cci.gainmap_sagged_pixel_overplotter
   :members:
