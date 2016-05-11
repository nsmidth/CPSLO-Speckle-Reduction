import sys, os
from classes_labeyrie import target,deconvolved
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np

binary = target()
reference = target()
deconv = deconvolved()

# Define filenames for each target
binary.fits.fileName = "/home/niels/Documents/FITS/KP336.fits"
reference.fits.fileName = "/home/niels/Documents/FITS/KP338.fits"

# Import each target
binary.fits.read(numDimensions=3,printInfo=False)
reference.fits.read(numDimensions=3,printInfo=False)

# Calculate PSD of each target
binary.psdCalc()
print("Binary PSD Calc Complete")
reference.psdCalc()
print("Reference PSD Calc Complete")

# Deconvolve reference and binary stars
deconv.psdDeconvolve(binary.psd.data,reference.psd.data,1e-12)

# Perform filtering on output star object
deconv.psdFilter(lpfRadius = 20, interference=False)

# Get autocorrelogram
deconv.acorrCalc()

plt.figure(figsize = (12,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(binary.fits.data[10], cmap=colormap)
plt.title("Binary Star Image")
plt.subplot(1,2,2)
plt.imshow(reference.fits.data[0], cmap=colormap)
plt.title("Reference Star Image")

plt.figure(figsize = (6,12), dpi = 100)
plt.subplot(2,1,1)
plt.imshow(deconv.psdFiltered.data, cmap=colormap)
plt.title("Filtered Deconvolved PSD")
plt.subplot(2,1,2)
plt.imshow(deconv.acorr.data, cmap=colormap)
plt.title("Autocorrelation")

plt.show()
