""" 

    Testing to see if the output of PS3's pre-processing is the same as 
	output of current preprocessing algorithm

    Place the preprocessed file in the same folder as the non-processed one,
   	with the suffix _PSD to indicate it is preprocessed

    Niels Smidth - 10/23/15
"""

# Import Modules
from speckle_fns import fits_import,preprocess

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fftshift

import tkinter as tk
from tkinter import filedialog

# Query user for FITS files
print("Place the preprocessed FITS file with suffix _PSD in")
print(" the same directory as the unprocessed FITS file")
print("Select a FITS file")
root = tk.Tk()  # Starting up Tkinter
root.withdraw()
filePath = filedialog.askopenfilename() # Filename dialog

filePathPSD = filePath.replace(".fits","_PSD.fits") # PS3 FITS File

# Open FITS Data
fitsData = fits_import(filePath)

# Pre-Process FITS Data
psdAvg = preprocess(fitsData)

# Opening PS3's Preprocessed FITS file
psdAvgPSD = fits_import(filePathPSD)

# Check if arrays are the same
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
    plt.imshow(np.log10(1 + fftshift(pyImage) ))
    plt.title('Python PSD')

    plt.subplot(1,2,2)
    plt.imshow(np.log10(1 + fftshift(ps3Image) ))
    plt.title('PS3 PSD')

    plt.show()
