# Calculates PSD of input FITS files

# Import Modules
import sys, os
from labeyrie_classes import target
import tkinter as tk
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt 

# Raw data target object
raw = target()

# Start up Tkinter
root = tk.Tk()
root.withdraw()

# Prompt for FITS file locations
fitsFileNames = filedialog.askopenfilenames(title="Select FITS files for processing")
# Debug
#fitsFileNames = []
#fitsFileNames.append('/home/niels/Desktop/KP330.fits')
#fitsFileNames.append('/home/niels/Desktop/KP331.fits')
#fitsFileNames.append('/home/niels/Desktop/KP332.fits')

# Loop through each fileName
for fitsFileName in fitsFileNames:

    # Import FITS file
    raw.fits.fileName = fitsFileName
    raw.fits.read(numDimensions=3,printInfo=False)

    # Create new filename
    psdFileName = os.path.splitext(fitsFileName)[0]
    psdFileName = psdFileName + "_PSD.fits"

    # If filename already exists, add "0" to the end of filename
    #  until no longer is a repeat
    # Safety condition for this loop : break out if more than 4 repeats
    fileRepeat = 0
    while( os.path.exists(psdFileName) ):
        # Check if too many repeated filenames
        if (fileRepeat > 3):
            print("Too many repeated filenames, delete some")
            break
        # Add a 0 to filename
        psdFileName = os.path.splitext(psdFileName)[0]
        psdFileName = psdFileName + "0.fits"
        fileRepeat += 1


    if (os.path.exists(psdFileName)): # if no good filename found
        # Don't generate file
        print("No processed file generated")


    # If good filename found
    else:
        print("Processing file: ", fitsFileName)
        print("Creating file: ", psdFileName)

        # Process FITS data
        raw.psdCalc()

        # Create new FITS file
        raw.psd.fileName = psdFileName
        raw.psd.write()

        #Print message for user
        print("Done processing ", psdFileName)
        print()
        print()

