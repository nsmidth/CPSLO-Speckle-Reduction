# Module for performing Labeyrie deconvolution of preprocessed 
#  FITS files

# Included modules
import sys
from labeyrieClasses import target,deconvolved, photometry, camsky
import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np
import math
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

# Write acorr to file
deconv.acorr.fileName = ''
deconv.acorr.fileName = filedialog.asksaveasfilename(defaultextension=".fits",
                                                     initialdir=(os.path.split(reference.psd.fileName+'/')))
if (deconv.acorr.fileName != ''):
    deconv.acorr.write()

# Displaying estimated observed and expected locations
# Calculate middle index of image
midpoint = np.shape(deconv.acorr.data)[0]/2
# Create object for observed and expected secondary locations
obs0 = camsky(midpoint=midpoint,delta=0,e=0.01166)
obs1 = camsky(midpoint=midpoint,delta=0,e=0.01166)
exp = camsky(midpoint=midpoint,delta=0,e=0.01166)

# Input expected secondary location to be displayed
## DEBUG################# KP336/338 Expected values
exp.sky.theta = 305
exp.sky.rho = 0.5
exp.sky2cam()
exp.cam.polar2cart()

# Testing centroid estimation
calcs.acorr = deconv.acorr.data # Moving acorr to photometry object
centroid = calcs.centroidEstimate()
# Save centroid locations to observed position objects
(obs0.cam.x,obs0.cam.y)=centroid[0]
(obs1.cam.x,obs1.cam.y)=centroid[1]

# Viewing estimated centroid locations
# Prepare the object for marking up with locations
calcs.acorr = deconv.acorr.data
calcs.acorrMarkedClear()
# Mark expected position
calcs.acorrMark(exp.cam.x.astype(np.uint16),exp.cam.y.astype(np.uint16),"o",(0,255,0))
# Mark observed position
calcs.acorrMark(obs0.cam.x.astype(np.uint16),obs0.cam.y.astype(np.uint16),"o",(255,0,0))
calcs.acorrMark(obs1.cam.x.astype(np.uint16),obs1.cam.y.astype(np.uint16),"o",(255,0,0))
# View Image
plt.figure()
plt.imshow(calcs.acorrMarked)
plt.show()