# Module for performing Labeyrie deconvolution of preprocessed 
#  FITS files

# Included modules
import sys
from labeyrieClasses import target,deconvolved, photometry, camsky
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
binary.psd.fileName = "/home/niels/Documents/FITS/KP336_PSD.fits"
# Import binary star PSD FITS
binary.psd.read()

# Prompt user for reference PSD file location
#reference.psd.fileName = filedialog.askopenfilename(title="Select REFERENCE FITS file")
reference.psd.fileName = "/home/niels/Documents/FITS/KP338_PSD.fits"
# Import reference star PSD FITS
reference.psd.read()

# Deconvolve reference and binary stars
deconv.psdDeconvolve(binary.psd.data,reference.psd.data,1e-12)

# Perform filtering on output star object
deconv.psdFilter(lpfRadius = 20, interference=True)

# Get autocorrelogram
deconv.acorrCalc()

# View Results
deconv.acorr.view(log=True)

# Estimating centroid locations
# Moving acorr to photometry object
calcs.acorr = deconv.acorr.data
calcs.acorrMarkedClear()

# Calculate middle index of image
midpoint = np.shape(calcs.acorr)[0]/2
# Create object for observed and expected secondary locations
obs = camsky(midpoint=midpoint,delta=0,e=0.01166)
exp = camsky(midpoint=midpoint,delta=0,e=0.01166)

# Input expected secondary location to be displayed
## DEBUG################# KP336/338 Expected values
exp.sky.theta = 305
exp.sky.rho = 0.5
exp.sky2cam()
exp.cam.polar2cart()

# Testing centroid estimation
centroid = calcs.centroidEstimate()
centroid = np.round(centroid).astype(np.uint16)

# Viewing estimated centroid locations
plt.figure()
# Mark expected position
calcs.acorrMark(exp.cam.x.astype(np.uint16),exp.cam.y.astype(np.uint16),"o",(0,255,0))
# Mark observed position
calcs.acorrMark(centroid[0,0],centroid[0,1],"+",(255,0,0))
calcs.acorrMark(centroid[1,0],centroid[1,1],"+",(255,0,0))
# View Image
plt.imshow(calcs.acorrMarked)
plt.show()