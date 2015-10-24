'''
    Function to easily create a simple FITS file

    Niels Smidth - 10/23/15
'''

import numpy as np
from astropy.io import fits


# filename :    Filename for FITS file to be generated, must be of
#               format "xxxx.fits"
# data :        Numpy Array to be used as FITS primary HDU data
def fits_generate(filename, data):
    
    # Create PrimaryHDU object with data
    hdu = fits.PrimaryHDU(data)

    # Create HDUList object w/ PrimaryHDU
    hdulist = fits.HDUList([hdu])

    # Write to new file
    hdulist.writeto(filename)
