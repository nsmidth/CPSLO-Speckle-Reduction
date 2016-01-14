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
dcnvlv = deconvolved()

# Prompt user for binary PSD file location
binary.psdFileName = filedialog.askopenfilename(title="Select BINARY FITS file")
# Import binary star PSD FITS
binary.psdImport()
#binary.psdView()

# Prompt user for reference PSD file location
reference.psdFileName = filedialog.askopenfilename(title="Select REFERENCE FITS file")
# Import reference star PSD FITS
reference.psdImport()
#reference.psdView()

# Deconvolve reference and binary stars
#dcnvlv.psdDeconvolve(binary.PSD,reference.PSD)

# Perform filtering on output star object
#dcnvlv.psdFilter(lpf = 15, interference = TRUE)

# Get autocorrelogram
#dcnvlv.acorrCalc()

# Return autocorrelogram as PNG or FITS
