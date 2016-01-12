# Module for performing Labeyrie deconvolution of preprocessed 
#  FITS files

# Included modules
import sys
from labeyrieClasses import target,deconvolved
import tkinter as tk
from tkinter import filedialog

# Start up Tkinter
root = tk.Tk()  
root.withdraw()

# Instantiate objects
binary = target()
reference = target()
dcnvlv = deconvolved()

# Prompt user for binary PSD file location
#binary.psdFileName = 
# Import binary star PSD FITS
#binary.psdImport()

# Prompt user for reference PSD file location
#reference.psdFileName = 
# Import reference star PSD FITS
#reference.psdImport()

# Deconvolve reference and binary stars
#dcnvlv.psdDeconvolve(binary.PSD,reference.PSD)

# Perform filtering on output star object
#dcnvlv.psdFilter(lpf = 15, interference = TRUE)

# Get autocorrelogram
#dcnvlv.acorrCalc()

# Return autocorrelogram as PNG or FITS
