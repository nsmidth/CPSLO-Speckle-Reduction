""" 

    Attempting to recreate PS3's Speckle Processing

    Niels Smidth - 10/16/15
"""

# Import Modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2,ifft2, fftshift

from astropy.io import fits

import tkinter as tk
from tkinter import filedialog

# Importing FITS Data
print("Please select a FITS file")
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()
HDUList = fits.open(file_path)

# Print FITS File Info & Headers
HDUList.info()
print()
print("Headers:")
print(repr(HDUList[0].header))

# Save data in FITS cube to local variable, then close FITS file
fitsData = HDUList[0].data
HDUList.close()

# Generate empty array the size of an image to be used to accumulate
#  PSD values before averaging. 
psdSum = np.zeros(fitsData.shape[1:3])

# Looping through all images in cube
for img in fitsData[:]:

    # FFT function requires little-endian data, so casting it
    imgStar = img.astype(float)

    ## Preprocessing of data
    # Take FFT of images, this gives us complex numbers
    fftStar = fft2(imgStar)

    # Calculate 2D power spectrum
    # This gives us only real values as 
    absStar = np.abs( fftStar)
    psdStar = absStar**2

    # Doing FFT shift on PSD, which moves low spatial frequencies to center
    psdStar = fftshift(psdStar)

    # Accumulate current PSD value
    psdSum += psdStar

# Divide by # of images to calculate average
psdAvg = psdSum/fitsData.shape[0]

# Do iFFT on PSD's, bringing back to spatial domain
# This should give us the autocorrelations of original images
acorrStar = ifft2(psdAvg)

# Taking iFFT of PSD (all real values) results in complex valued output
#  Must view the magnitude of the output
#  Doing FFTshift to move eyes of autocorrelation near center
acorrStar = np.abs(fftshift(acorrStar))
# Taking iFFT of PSD (all real values) results in complex valued output
#  Must view the magnitude of the output

# View images

plt.figure(num=1,figsize=(9,3),dpi=120)
plt.subplot(1,3,1)
plt.imshow(imgStar)
plt.title('Ex Star Image')

plt.subplot(1,3,2)
plt.imshow(np.log10( psdAvg ))
plt.title('Avg Star PSD')

plt.subplot(1,3,3)
plt.imshow(np.abs(acorrStar))
plt.title('Autocorrelation')

plt.show()


# Debug End of Code
print("Done Running")

