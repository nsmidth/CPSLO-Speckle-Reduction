""" This is a test of importing an artificial speckle image
    and performing some frequency domain analysis on it

    Niels Smidth - 10/16/15
"""

# Import Modules
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from scipy.fftpack import fft2,ifft2,fftshift

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

# Calculate difference of log10(PSD's), this is equivalent
#  to deconvolving the two
psdDeconvolved = np.log10(psdDoubleStar) - np.log10(psdSingleStar)
psdDeconvolved = 10**psdDeconvolved

# Doing FFT shift on PSD, which moves low spatial frequencies to center
psdSingleStar = fftshift(psdSingleStar)
psdDoubleStar = fftshift(psdDoubleStar)
psdDeconvolved = fftshift(psdDeconvolved)

# Do iFFT on PSD's, bringing back to spatial domain
# This should give us the autocorrelations of original images
acorrSingleStar = ifft2(psdSingleStar)
acorrDoubleStar = ifft2(psdDoubleStar)
acorrDeconvolved = ifft2(psdDeconvolved)

# Taking iFFT of PSD (all real values) results in complex valued output
#  Must view the magnitude of the output
#  Doing FFTshift to move eyes of autocorrelation near center
acorrSingleStar = np.abs(fftshift(acorrSingleStar))
acorrDoubleStar = np.abs(fftshift(acorrDoubleStar))
acorrDeconvolved = np.abs(fftshift(acorrDeconvolved))

# View images
## Figure 1 : Single Star
plt.figure(num=1, figsize=(9,3), dpi = 120)
plt.subplot(1,3,1)
plt.imshow(imgSingleStar, cmap=plt.cm.Greys)
plt.title('Single Star Image')

plt.subplot(1,3,2)
plt.imshow(np.log10( psdSingleStar ), cmap=plt.cm.Greys)
plt.title('Single Star PSD')

plt.subplot(1,3,3)
plt.imshow( acorrSingleStar , cmap=plt.cm.Greys)
plt.title('Single Star Autocorrelation')

## Figure 2 : Double Star
plt.figure(num=2, figsize=(9,3), dpi = 120)
plt.subplot(1,3,1)
plt.imshow(imgDoubleStar, cmap=plt.cm.Greys)
plt.title('Double Star Image')

plt.subplot(1,3,2)
plt.imshow(np.log10( psdDoubleStar ), cmap=plt.cm.Greys)
plt.title('Double Star PSD')

plt.subplot(1,3,3)
plt.imshow( acorrDoubleStar , cmap=plt.cm.Greys)
plt.title('Double Star Autocorrelation')

## Figure 3 : Deconvolution
plt.figure(num=3,figsize=(6,3),dpi=120)
plt.subplot(1,2,1)
plt.imshow(np.log10(psdDeconvolved), cmap=plt.cm.Greys)
plt.title('Deconvolution PSD')

plt.subplot(1,2,2)
plt.imshow(acorrDeconvolved, cmap=plt.cm.Greys)
plt.title('Deconvolution Autocorrelation')

plt.show()

# Debug End of Code
print("Done Running")
