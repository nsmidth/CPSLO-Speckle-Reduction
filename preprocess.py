# Calculates PSD of input FITS files

# Import Modules
import sys, os
from labeyrieClasses import target
import tkinter as tk
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt 

# Raw data target object
raw = target()

# Prompt for FITS file locations
fitsFileNames = filedialog.askopenfilenames(title="Select FITS files for processing")


# Loop through each fileName
for fitsFileName in fitsFileNames:
    # Checking if fitsFileName is a .fits
    if (os.path.splitext(fitsFileName)[1] == ".fits"): # If filetype is FITS, preprocess

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
            #print("Processing file: ", fitsFileName)
            #print("Creating file: ", psdFileName)
                      
            # Import FITS file
            raw.fitsFileName = fitsFileName
            raw.psdFileName = psdFileName
            raw.fits = raw.fitsImport()
            raw.fitsView(1)

            # Preprocess FITS data
            #raw.preprocess()
            
            # Create new FITS file 
            #raw.fitsExport()
            
            #Print message for user
            #print("Done processing ", psdFileName)

    else: # If filetype is not FITS, don't preprocess
        # Print message for user
        print("The following file is not .fits: ", fitsFileName)

print()
print()
    

    

#Debug

