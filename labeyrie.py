# Module for performing Labeyrie deconvolution of preprocessed 
#  FITS files

# Included modules
import sys
from labeyrieClasses import target,deconvolved
import tkinter as tk
from tkinter import filedialog
import os, sys

# Start up Tkinter
root = tk.Tk()
root.withdraw()

# Instantiate objects
binary = target()
reference = target()
soln = deconvolved()

# Prompt user for binary PSD file location
#binary.psdFileName = filedialog.askopenfilename(title="Select BINARY FITS file")
binary.psdFileName = "/home/niels/Desktop/KP330_PSD.fits" #DEBUG
# Import binary star PSD FITS
binary.psdImport()
#binary.psdView()

# Prompt user for reference PSD file location
#reference.psdFileName = filedialog.askopenfilename(title="Select REFERENCE FITS file")
reference.psdFileName = "/home/niels/Desktop/KP331_PSD.fits" #DEBUG
# Import reference star PSD FITS
reference.psdImport()
#reference.psdView()

# Deconvolve reference and binary stars
soln.psdDeconvolve(binary.psd,reference.psd,1e-12)
soln.psdView()

# Perform filtering on output star object
#soln.psdFilter(lpf = 15, interference = TRUE)
soln.psdFilter(lpfRadius = 25, interference=True)
soln.psdFilteredView()

# Get autocorrelogram
#soln.acorrCalc()

# Return autocorrelogram as PNG or FITS
