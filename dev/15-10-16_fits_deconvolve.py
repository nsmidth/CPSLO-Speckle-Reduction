"""
    Testing deconvolution of Preprocessed FITS files

    Niels Smidth - 10/16/15
"""

# Import Modules
import matplotlib.pyplot as plt
from scipy.fftpack import fftshift
import numpy as np
from speckle_fns import postprocess,fits_import
from speckle_fns import deconv0,deconv1
import tkinter as tk
from tkinter import filedialog
import os as os
import sys

# Start up Tkinter
root = tk.Tk()   
root.withdraw()

# Query user for FITS files to be deconvolved
#filePathDouble = filedialog.askopenfilename(title="Please select the preprocessed Double Star FITS file")
#filePathSingle = filedialog.askopenfilename(title="Please select the preprocessed Single Star FITS file")

# Debug : Hardcoded filenames
filePathDouble = "/home/niels/Desktop/KP330_PSD.fits"
filePathSingle = "/home/niels/Desktop/KP331_PSD.fits"

# Confirm files are .fits files
#if( os.path.split(filePathDouble)[1] != ".fits" ):
#    print(

# Import FITS file data
fitsDoubleStar = fits_import(filePathDouble)
fitsSingleStar = fits_import(filePathSingle)

# Confirm that inputted data is a 2 dimensional array
#if ( len(np.shape(fitsDoubleStar)) != 2 ):

# Convert FITS data to floats
psdDoubleStar = fitsDoubleStar.astype(float)
psdSingleStar = fitsSingleStar.astype(float)

# Perform deconvolution on images
constant = 1E-15
psdDeconv1 = deconv1(psdDoubleStar, psdSingleStar, constant)

# Perform postprocessing on images
acorrSingleStar = postprocess( psdSingleStar )
acorrDoubleStar = postprocess( psdDoubleStar )
acorrDeconv1 = postprocess( psdDeconv1 )

# View images
printImg = input("Print star images (y/n)?   ")
if (printImg.lower() == "y"):
    # Figure 1 shows PSD plots
    plt.figure(num=1, figsize=(12,3), dpi=128)
    plt.subplot(1,3,1)
    plt.title("PSD Reference")
    plt.imshow(np.log10(fftshift(psdSingleStar)))
    plt.colorbar()

    plt.subplot(1,3,2)
    plt.title("PSD Binary")
    plt.imshow(np.log10(fftshift(psdDoubleStar)))
    plt.colorbar()

    plt.subplot(1,3,3)
    plt.title("PSD Deconv")
    plt.imshow(np.log10(fftshift(psdDeconv1)))
    plt.colorbar()

    # Figure 2 shows Acorr Plots
    plt.figure(num=2, figsize=(12,3), dpi=128)
    plt.subplot(1,3,1)
    plt.title("Acorr Reference")
    plt.imshow(acorrSingleStar)
    plt.colorbar()

    plt.subplot(1,3,2)
    plt.title("Acorr Binary")
    plt.imshow(acorrDoubleStar)
    plt.colorbar()

    plt.subplot(1,3,3)
    plt.title("Acorr Deconv")
    plt.imshow(acorrDeconv1)
    plt.colorbar()


    plt.show(block=False)

# Debug End of Code
print("Done Running")
