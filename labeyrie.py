# Module for performing Labeyrie deconvolution of preprocessed 
#  FITS files

# Included modules
import sys
from labeyrieClasses import target,deconvolved, photometry
import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np
from tkinter import filedialog
import os, sys

# Start up Tkinter
root = tk.Tk()
root.withdraw()

# Instantiate objects
binary = target()
reference = target()
deconv = deconvolved()
calcs = photometry()

# Prompt user for binary PSD file location
#binary.psd.fileName = filedialog.askopenfilename(title="Select BINARY FITS file")
binary.psd.fileName = "/home/niels/Desktop/KP330_PSD.fits"
# Import binary star PSD FITS
binary.psd.read()

# Prompt user for reference PSD file location
#reference.psd.fileName = filedialog.askopenfilename(title="Select REFERENCE FITS file")
reference.psd.fileName = "/home/niels/Desktop/KP331_PSD.fits"
# Import reference star PSD FITS
reference.psd.read()

# Deconvolve reference and binary stars
deconv.psdDeconvolve(binary.psd.data,reference.psd.data,1e-12)

# Perform filtering on output star object
deconv.psdFilter(lpfRadius = 25, interference=True)

# Get autocorrelogram
deconv.acorrCalc()

# View Results
#deconv.acorr.view()

# Estimating centroid locations
# Moving acorr to photometry object
calcs.acorr = deconv.acorr.data
calcs.acorrMarkedClear()

# Testing centroid estimation
centroid = calcs.centroidEstimate()
centroid = np.round(centroid).astype(np.uint16)
print(centroid)

# Viewing estimated centroid locations
plt.figure()
calcs.acorrMark(centroid[0,0],centroid[0,1],"+",(255,0,0))
calcs.acorrMark(centroid[1,0],centroid[1,1],"+",(255,0,0))
plt.imshow(calcs.acorrMarked)
plt.show()