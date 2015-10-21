""" This is a test of importing an artificial speckle image
    and performing some frequency domain analysis on it

    Niels Smidth - 10/16/15
"""

# Import Modules
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from scipy.fftpack import fft2,ifft2

# Import images
imgSingleStar = mpimg.imread('10-15 Single Speckle Test.png')
imgSingleStar = imgSingleStar[:,:,1] # Only care about one channel
imgDoubleStar = mpimg.imread('10-15 Binary Speckle Test.png')
imgDoubleStar = imgDoubleStar[:,:,1] # Only care about one channel

## Preprocessing of data
# Take FFT of images, this gives us complex numbers
fftSingleStar = fft2(imgSingleStar)
fftDoubleStar = fft2(imgDoubleStar)

# Calculate 2D power spectrum
# This gives us only real values as 
absSingleStar = np.abs( fftSingleStar)
absDoubleStar = np.abs( fftDoubleStar)

psdSingleStar = absSingleStar**2
psdDoubleStar = absDoubleStar**2

# For a FITS cube, would calculate averages of all the PSDs

## Postprocessing of data
# Calculate quotient of double star and reference single star (deconvolution)
#deconvolved = psdDoubleStar/psdSingleStar

# Calculate ifft of deconvolved data = autocorrelation
#autocorr = ifft2(deconvolved)

# View images
## Figure 1 : Single Star
plt.figure(1)
plt.subplot(2,1,1)
plt.imshow(imgSingleStar, cmap=plt.cm.Greys)
plt.title('Single Star Image')


plt.subplot(2,1,2)
plt.imshow(np.log10( psdSingleStar ), cmap=plt.cm.Greys)
plt.title('Single Star PSD')


## Figure 2 : Double Star
plt.figure(2)
plt.subplot(2,1,1)
plt.imshow(imgDoubleStar, cmap=plt.cm.Greys)
plt.title('Double Star Image')


plt.subplot(2,1,2)
plt.imshow(np.log10( psdDoubleStar ), cmap=plt.cm.Greys)
plt.title('Double Star PSD')


plt.show()

# Debug End of Code
print("Done Running")
