""" 

    Attempting to recreate PS3's Speckle Processing

    Niels Smidth - 10/16/15
"""

# Import Modules
from speckle_fns import fits_import,preprocess,postprocess,print_star

import tkinter as tk
from tkinter import filedialog

# Query user for FITS file name
print("Please select a FITS file")
root = tk.Tk()  # Starting up Tkinter
root.withdraw()
file_path = filedialog.askopenfilename() # Filename dialog

# Import FITS data cube
fitsData = fits_import(file_path)

# Perform pre-processing on FITS data
psdAvg = preprocess(fitsData)

# Perform post-processing on pre-processed data
acorrStar = postprocess(psdAvg)

# Print star if user desires
printImg = input("Print star images (y/n)?   ")
if (printImg.lower() == "y"):
    print_star(fitsData[0], psdAvg, acorrStar)

# Debug End of Code
print("Done Running")

