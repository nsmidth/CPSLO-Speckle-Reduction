""" This is a test of importing a speckle image
    and performing some frequency domain analysis on it

    Niels Smidth - 10/16/15
"""

# Import Modules
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from scipy.fftpack import fft2,ifft2, fftshift
import tkinter as tk
from tkinter import filedialog

# Import images
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()
imgStar = mpimg.imread(file_path)
imgStar = imgStar[:,:,1] # Only care about one channel

## Preprocessing of data
# Take FFT of images, this gives us complex numbers
fftStar = fft2(imgStar)

# Calculate 2D power spectrum
# This gives us only real values as 
absStar = np.abs( fftStar)
psdStar = absStar**2

# Doing FFT shift on PSD, which moves low spatial frequencies to center
psdStar = fftshift(psdStar)

# Do iFFT on PSD's, bringing back to spatial domain
# This should give us the autocorrelations of original images
acorrStar = ifft2(psdStar)

# Taking iFFT of PSD (all real values) results in complex valued output
#  Must view the magnitude of the output
#  Doing FFTshift to move eyes of autocorrelation near center
acorrStar = np.abs(fftshift(acorrStar))
# Taking iFFT of PSD (all real values) results in complex valued output
#  Must view the magnitude of the output

# View images
plt.figure(num=1,figsize=(9,3),dpi=120)
plt.subplot(1,3,1)
plt.imshow(imgStar, cmap=plt.cm.Greys)
plt.title('Star Image')


plt.subplot(1,3,2)
plt.imshow(np.log10( psdStar ), cmap=plt.cm.Greys)
plt.title('Star PSD')

plt.subplot(1,3,3)
plt.imshow( np.abs(acorrStar) , cmap=plt.cm.Greys)
plt.title('Star Autocorrelation')

plt.show()

# Debug End of Code
print("Done Running")
