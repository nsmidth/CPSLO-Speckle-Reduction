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
binary.psd.fileName = filedialog.askopenfilename(title="Select BINARY FITS file")
# Import binary star PSD FITS
binary.psd.read()
binary.psd.view(log=True)

# Prompt user for reference PSD file location
reference.psd.fileName = filedialog.askopenfilename(title="Select REFERENCE FITS file")
# Import reference star PSD FITS
reference.psd.read()
reference.psd.view(log=True)

# Deconvolve reference and binary stars
soln.psdDeconvolve(binary.psd.data,reference.psd.data,1e-12)

# Perform filtering on output star object
soln.psdFilter(lpfRadius = 25, interference=True)

# Get autocorrelogram
soln.acorrCalc()

# View Results
soln.psdFiltered.view(log=True)
soln.acorr.view()
