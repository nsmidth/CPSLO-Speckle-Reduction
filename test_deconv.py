#
#
# Testing several different deconvolution methods
# Viewing the transfer function for each deconvolution method
#
#

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft2,ifft2,fftshift
from speckle_fns import deconv0,deconv1,fits_import

# Import the reference star PSD that would be used in deconvolution
refStar = fits_import("/home/niels/Desktop/KP331_PSD.fits")

# Create an array of 1's to be used with deconvolution functions
#  as the double star array
unityStar = np.zeros( (512,512) )
unityStar.fill(1)

# Show reference star PSD
plt.figure()
plt.imshow(np.log10(fftshift(refStar)))
plt.colorbar()
#plt.set_cmap('copper')
plt.title("Reference Star PSD")

# Test deconv0 function
print("d0 = deconv0 function")
d0 = deconv0(unityStar, refStar, refStar.min()*0.01)

plt.figure()
plt.imshow(np.log10(fftshift(d0)))
plt.colorbar()
plt.title("d0")

# Test deconv1 function
print("d1 = deconv1 function")
d1 = deconv1(unityStar, refStar, refStar.min()*1E-2)

plt.figure()
plt.imshow(np.log10(fftshift(d1)))
plt.colorbar()
plt.title("d1")

# Test "pseudo-inverse" from Mathworks Image Deblurring Intro
print("d2 = pseudo-inverse from Mathworks Image Deblurring Intro")
d2 = refStar/(refStar**2 + (1E-2)**2)

plt.figure()
plt.imshow(np.log10(fftshift(d2)))
plt.colorbar()
plt.title("d2")

# Test wiener filter autocorrelation model approximation
print("d3 = wiener filter autoccorrelation model approximation")
# Importing all images to get image/noise statistics
fitsData = fits_import("/home/niels/Desktop/KP331.fits")
fitsAvg = np.average(fitsData, axis=0)

# Looking at only the noise of image


# Show images
plt.show(block=False)
