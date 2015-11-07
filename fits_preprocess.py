"""

    Preprocess FITS files.
    User selects several FITS files to be procesed.
    Average PSD is calculated for each file.
    Average PSD is saved into a new FITS file with _PSD suffix

    Niels Smidth - 11/6/15

"""

# Import Modules
from speckle_fns import fits_import,preprocess,fits_generate
import numpy as np
import os
import time
import tkinter as tk
from tkinter import filedialog

# Start up Tkinter
root = tk.Tk()  
root.withdraw()

# Query user for FITS files to be preprocessed
file_path = filedialog.askopenfilenames(title="Please select the FITS files to be preprocessed")

# Loop through each file
for file in file_path:
    # Checking if file is a .fits
    fileExt = file.split('.')   # Splitting by periods to find file extension
    indexExt = ( np.shape(fileExt)[0] - 1 ) # File ext is last string
    fileExt = fileExt[indexExt]

    if (fileExt == "fits"): # If filetype is FITS, preprocess
        # Print message for user
        print("Processing file: ",file)

        # Create new filename      
        filePSD = os.path.splitext(file)[0]
        filePSD = filePSD + "_PSD.fits"
        
        # If filename already exists, add "0" to the end of filename
        #  until no longer is a repeat
        # Safety condition for this loop : break out if more than 4 repeats
        fileRepeat = 0
        while( os.path.exists(filePSD) ):
            # Check if too many repeated filenames
            if (fileRepeat > 3):
                print("Too many repeated filenames, delete some")
                break
            # Add a 0 to filename
            filePSD = os.path.splitext(filePSD)[0]
            filePSD = filePSD + "0.fits"
            fileRepeat += 1


        if (os.path.exists(filePSD)): # if no good filename found
            # Don't generate file 
            print("No processed file generated")
            
        else: # if good filename found
            # Import FITS file
            fitsData = fits_import(file)

            # Preprocess FITS data
            fitsPSD = preprocess(fitsData)
            
            # Create new FITS file 
            fits_generate(filePSD, fitsPSD)
            
            #Print message for user
            print("Created ", filePSD)

        print()
        print()

    else: # If filetype is not FITS, don't preprocess
        # Print message for user
        print("The following file is not .fits: ", file)
        print()
        print()
    

    

    


