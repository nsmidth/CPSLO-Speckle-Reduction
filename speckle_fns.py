""" 

    Attempting to recreate PS3's Speckle Processing

    Niels Smidth - 10/16/15
"""

# Import Modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2,ifft2, fftshift

from astropy.io import fits

# fits_import()
# Queries user for FITS file to import, imports FITS file,
#  optionally prints info/headers, and returns data. Assumes
#  that data is contained in primary HDU
#
# Args:
#  file_path : filepath of FITS file to be opened
#  printInfo : set to TRUE to print FITS file info
#  printHeaders : set to TRUE to print FITS file headers
#
# Returns data block from primary HDU
#
def fits_import(file_path, printInfo = True, printHeaders = False):
    
    HDUList = fits.open(file_path)  # Open FITS Data

    # Print FITS File Info & Headers
    if (printInfo == True):
        HDUList.info()
        print()

    if (printHeaders == True):
        print("Headers:")
        print(repr(HDUList[0].header))

    # Save data in FITS cube to local variable, then close FITS file
    fitsData = HDUList[0].data
    HDUList.close()

    # Return FITS data
    return fitsData

# preprocess()
# Performs a similar pre-processing on FITS images as PS3
#  Each Image -> FFT (spectrum) -> magnitude^2 (PSD) -> Averaged
#
# Args:
#  fitsData : Array of data to be pre-processed. Assume it is in
#   following form, 100x512x512 for 100 square images
#
# Returns averaged PSD array
def preprocess(fitsData):
    # Checking if FITS data is an array of images
    if (len(fitsData.shape) == 3):
        # Generate empty array the size of an image to be used to accumulate
        #  PSD values before averaging. 
        psdSum = np.zeros(fitsData.shape[1:3])

        # Looping through all images in cube
        for index,img in enumerate(fitsData[:]):

            # Print current file being processed
            print("Processing Image #: ",(index+1))

            # FFT function requires little-endian data, so casting it
            imgStar = img.astype(float)

            ## Preprocessing of data
            # Take FFT of images, this gives us complex numbers
            fftStar = fft2(imgStar)

            # Calculate 2D power spectrum
            # This gives us only real values as 
            absStar = np.abs( fftStar)
            psdStar = absStar**2

            # Accumulate current PSD value
            psdSum = np.add(psdSum,psdStar)

        # Divide by # of images to calculate average
        psdAvg = np.divide(psdSum,fitsData.shape[0])

    #Otherwise if FITS data is only one image
    elif (len(fitsData.shape) == 2):
        # FFT function requires little-endian data, so casting it
        imgStar = fitsData.astype(float)

        ## Preprocessing of data
        # Take FFT of images, this gives us complex numbers
        fftStar = fft2(imgStar)

        # Calculate 2D power spectrum
        # This gives us only real values as 
        absStar = np.abs( fftStar)
        psdAvg = absStar**2

    return psdAvg

# postprocess()
#
# Brings the averaged PSD image back to an autocorrelation image
#
# Args:
#  psdAvg : PSD to be converted
#
# Returns autocorrelation image
def postprocess( psdAvg ):
    # Do iFFT on PSD's, bringing back to spatial domain
    # This should give us the autocorrelations of original images
    acorrStar = ifft2(psdAvg)

    # Taking iFFT of PSD (all real values) results in complex valued output
    #  Must view the magnitude of the output
    #  Doing FFTshift to move eyes of autocorrelation near center
    acorrStar = np.abs(fftshift(acorrStar))
    # Taking iFFT of PSD (all real values) results in complex valued output
    #  Must view the magnitude of the output

    return acorrStar

# print_star()
#
# Display image of all the important data for a star
#
# Args:
#  imgStar : A spatial domain image of star for reference
#  psdAvg : Averaged PSD of star
#  acorrStar : Autocorrelation of a star
def print_star(imgStar, psdAvg, acorrStar):
    # View images
    plt.figure(figsize=(9,3),dpi=120)
    plt.subplot(1,3,1)
    plt.imshow(imgStar)
    plt.title('Ex Star Image')

    plt.subplot(1,3,2)
    plt.imshow(fftshift( np.log10(1+psdAvg)) )
    plt.title('Avg Star PSD')


    plt.subplot(1,3,3)
    plt.imshow(np.abs(acorrStar))
    plt.title('Autocorrelation')
    plt.colorbar()

    plt.show()

# fits_generate()
#
# Args:
#  filename : Filename for FITS file to be generated, must be of
#              format "xxxx.fits"
#  data : Numpy Array to be used as FITS primary HDU data
def fits_generate(filename, data):
    
    # Create PrimaryHDU object with data
    hdu = fits.PrimaryHDU(data)

    # Create HDUList object w/ PrimaryHDU
    hdulist = fits.HDUList([hdu])

    # Write to new file
    hdulist.writeto(filename)

# deconv0()
#
# Deconvolve a single reference star from binary star.
#  This method replaces matching zeros in single/double 
#  star PSDs with small value (to avoid divide by zero).
#  Then the double star PSD divided by single star PSD
#
# Args:
#  psdDoubleStar : numpy array of average PSD of double star
#  psdSingleStar : numpy array of average PSD of single star
#
# Returns deconvolved PSD
def deconv0(psdDoubleStar, psdSingleStar, zeroSub = 0.001):

    # In double star PSD, replace pixels that are 0 in both double and
    #  single star PSDs
#    psdDoubleStar[ (psdSingleStar == 0) and (psdDoubleStar == 0) ] = zeroSub

    # Replace all 0's in single star PSD
    psdSingleStar[psdSingleStar == 0] = zeroSub

    # Divide double PSD by reference single PSD to deconvolve
    psdDeconv = np.divide(psdDoubleStar, psdSingleStar)

    return psdDeconv
