""" 

    Testing to see if the output of PS3's pre-processing is the same as 
	output of current preprocessing algorithm

    Place the preprocessed file in the same folder as the non-processed one,
   	with the suffix _PSD to indicate it is preprocessed

    Niels Smidth - 10/23/15
"""

# Import Modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2,ifft2, fftshift

from astropy.io import fits

import tkinter as tk
from tkinter import filedialog

# Importing FITS Data
print("Place the preprocessed FITS file with suffix _PSD in")
print(" the same directory as the unprocessed FITS file")
print("Select a FITS file")
root = tk.Tk()  # Starting up Tkinter
root.withdraw()
# Filename dialog
filePath = filedialog.askopenfilename()
# Filename for preprocessed file
filePathPSD = filePath.replace(".fits","_PSD.fits") 


## Preprocessing FITS file
# Open FITS Data
HDUList = fits.open(filePath)  

# Print FITS File Info
HDUList.info()

# Save data in FITS cube to local variable, then close FITS file
fitsData = HDUList[0].data
HDUList.close()

# Generate empty array the size of an image to be used to accumulate
#  PSD values before averaging. 
psdSum = np.zeros(fitsData.shape[1:3])

# Looping through all images in cube
for index,img in enumerate(fitsData[:]):

    # Print current file being processed
    print("Processing Image #: ",(index+1))

    # FFT function requires little-endian data, so casting it
    imgStar = img.astype(float)

    ## Preprocessing of data
    # Take FFT of images, this gives us complex numbers
    fftStar = fft2(imgStar)

    # Calculate 2D power spectrum
    # This gives us only real values as 
    absStar = np.abs( fftStar)
    psdStar = absStar**2

    # Accumulate current PSD value
    psdSum = np.add(psdSum,psdStar)

# Divide by # of images to calculate average
psdAvg = np.divide(psdSum,fitsData.shape[0])

## Opening PS3's Preprocessed FITS file
# Open FITS Data
HDUList = fits.open(filePathPSD)

# Print FITS Info
HDUList.info()

# Save data in FITS cube to local variable, then close FITS file
psdAvgPSD = HDUList[0].data
HDUList.close()

if (np.array_equal(psdAvg, psdAvgPSD)):
    print("Preprocessed Files are SAME")
else:
    print("Preprocessed Files are NOT SAME")

# Debug End of Code
print("Done Running")

def displayImages (pyImage, ps3Image):
    # View images
    plt.figure(num=1,figsize=(9,3),dpi=120)
    plt.subplot(1,2,1)
    plt.imshow(np.log10( fftshift(pyImage) ))
    plt.title('Python PSD')

    plt.subplot(1,2,2)
    plt.imshow(np.log10( fftshift(ps3Image) ))
    plt.title('PS3 PSD')

    plt.show()
