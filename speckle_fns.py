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

        # Normalizing FFT
        # Divide by (#pixels)^2 after calculating PSD rather than
        #  each FFT(image) by #pixels in order to execute quicker
        psdAvg = np.divide(psdAvg, (psdAvg.size)**2)

    #Otherwise if FITS data is only one image
    elif (len(fitsData.shape) == 2):
        # FFT function requires little-endian data, so casting it
        imgStar = fitsData.astype(float)

        ## Preprocessing of data
        # Take FFT of images, this gives us complex numbers
        fftStar = fft2(imgStar)

        # Normalizing FFT
        fftStar = np.divide(fftStar,fftStar.size)

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
    # Because we normalized after FFT by multiplying by 1/N^2, and ifft
    #  function does this as well, we need to multiply by N^2 before ifft
    #  to prevent performing normalization twice
    psdAvg = psdAvg*(psdAvg.size)
    
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
#  This method replaces zeros in single star PSD and 
#  corresponding zeros in double star PSD with small
#  value (to avoid divide by zero). Then the double
#  star PSD divided by single star PSD
#
# Args:
#  psdDoubleStar : numpy array of average PSD of double star
#  psdSingleStar : numpy array of average PSD of single star
#
# Returns deconvolved PSD
def deconv0(psdDoubleStar, psdSingleStar, zeroSub = 0.001):

    # Modified Single and Double Star arrays (zeros removed)
    psdSingleStarMod = np.array(psdSingleStar)
    psdDoubleStarMod = np.array(psdDoubleStar)

    # Find values of 0 in psdSingleStar
    zeroIndex = np.where( psdSingleStarMod == 0 )
    # Set all 0's to zeroSub
    psdSingleStarMod[zeroIndex] = zeroSub

    # If the location corresponding to a 0 in psdSingleStar also has a 0 in
    #  psdDoubleStar, also change that to zeroSub
    for [index, value] in enumerate(psdDoubleStarMod[zeroIndex]):
        if ( psdDoubleStarMod[ zeroIndex[0][index], zeroIndex[1][index] ] == 0 ):
            psdDoubleStarMod[ zeroIndex[0][index], zeroIndex[1][index] ] = zeroSub

    # Divide double PSD by modified reference single PSD to deconvolve
    psdDeconv = np.divide(psdDoubleStarMod, psdSingleStarMod)    

    return psdDeconv

# deconv1()
#
# Deconvolve a single reference star from binary star.
#  This method adds a small value to single star image to
#  prevent divide by zero problems
#
# Args:
#  psdDoubleStar : numpy array of average PSD of double star
#  psdSingleStar : numpy array of average PSD of single star
#
# Returns deconvolved PSD
def deconv1(psdDoubleStar, psdSingleStar, constant = 0.001):

    # Divide double PSD by modified reference single PSD to deconvolve
    psdDeconv = np.divide(psdDoubleStar, (psdSingleStar+constant))    

    return psdDeconv

# circ_filter1()
#
# Create a frequency domain filter image
#  Rotates a 1D filter function around the (0,0) point
#  in order to turn it into a 2D filter
#
# Args:
#  radius : Radius of filter, all zeros outside this radius
#  size : Dim of square filter, set as same as image to be filtered
#  filter_fn : Takes single argument (radius), returns a single
#   number, the freq response at that radius/freq
#
# Returns filter image
def circ_filter1(radius, size, filter_fn):
    
    # Returns radius of pixel, as given in x and y coordinates
    def pixel_radius(x,y):
        radius = np.power( (np.power(x,2)+np.power(y,2)), 0.5 )
        return radius
    
    # Check if filter's radius is too large for defined size
    if (radius > (size/2.0)):
        print("Radius is larger than half the image size")
        sys.exit()
    
    # Center of FFT image indices, where filter will be centered.
    # Ex: For a 512,512 image, the circle should be centered at 256,256
    center = (size)/2.0
    
    # Image indices for looping through area of image with 
    #  the filter image
    edges = {"start":0,
             "end":0}
    edges["start"]=np.floor(center-radius)
    edges["end"]=np.ceil(center+radius)
    
    # Create Empty 2D array to hold pixel values
    image = np.zeros((size,size))
    
    # Loop through all pixels within square to draw circle
    # Loop through all rows
    for row in np.arange(edges['start'],edges['end']+1):
        # Loop through all columns of row
        for column in np.arange(edges['start'],edges['end']+1):
            # Check if pixel is within circle's radius
            pixelDist = pixel_radius( (row-center),(column-center) )
            if (pixelDist < radius):
                image[row,column] = filter_fn(pixelDist)
    
    return image


# filter_interference()
#
# Perform interference filtering on image. Center pixels (512/2=256) are
#  replaced with average of two adjecant pixels
#
# Args:
#  imgF : frequency domain image to be filtered
#
# Returns filtered image
#
def filter_interference(imgF):
   
    # Copy image to be filtered
    imgF_int = np.array(imgF)
    
    # Critical dimensions of image
    image_size = imgF_int.shape[0]
    int_index = image_size/2.0
    bad_pixels = np.array((int_index + 1, int_index - 1))

    # Perform filtering on vertical interference
    for row in np.arange(image_size):
        good_pixels_avg = np.average(imgF[row,(bad_pixels[0], bad_pixels[1])])
        imgF_int[row,256]=good_pixels_avg       
    
    # Perform filtering on horizontal interference
    for column in np.arange(512):
        good_pixels_avg = np.average(imgF[(bad_pixels[0], bad_pixels[1]),column])
        imgF_int[256,column]=good_pixels_avg  
        
    return imgF_int

# filter_lpf()
#
# Performs a gaussian lowpass filter on an input image, of selectable radius
#
# Args:
#  imgF : Frequency domain image to be filtered
#  lpf_rad : std_dev of gaussian =~ radius of filter
#
# Returns filtered image
#
def filter_lpf(imgF, lpf_rad):

    image_size = imgF.shape[0]
    circle_radius = (image_size/2)-2
    
    # Define gaussian function
    def filter_fn_lpf(radius):
        return np.exp(-(np.power(radius,2)/(2*np.power(lpf_rad,2))))

    # Calculate 2D Gaussian filter
    lpfF = circ_filter1(radius = circle_radius, size = image_size, filter_fn = filter_fn_lpf)
    
    # Return filtered fn
    return np.multiply(lpfF,imgF)

# filter_hpf()
#
# Performs a gaussian highpass filter on an input image, of selectable radius & amplitude
#
# Args:
#  imgF : Frequency domain image to be filtered
#  hpf_rad : std_dev of gaussian =~ radius of filter
#  hpf_amplitude : Number between 0 and 1. Limits minimum of filter.
#   A value of 1 gives a minimum of 0 in the LPF
#   A value of 0.8 gives a minimum greater than 0 in the HPF 
#
# Returns filtered image
## Perform HPF
def filter_hpf(imgF, hpf_rad, hpf_amplitude):

    image_size = imgF.shape[0]
    circle_radius = (image_size/2)-2
    
    # Define LP gaussian function
    def filter_fn_hpf(radius):  
        return np.multiply(hpf_amplitude,np.exp(-(np.power(radius,2)/(2*np.power(hpf_rad,2)))))

    # Calculate 2D Gaussian filter
    # Generate lowpass filter to serve as inverse of highpass filter
    hpf = circ_filter1(radius = circle_radius, size = image_size, filter_fn = filter_fn_hpf)

    # Generate array of 1's to be subtracted from
    ones = np.zeros((image_size,image_size))
    ones.fill(1)

    # Calculate 1 - LP, which is HP
    hpfF = np.subtract(ones,hpf)
    
    # Return Filtered
    return np.multiply(hpfF,imgF)
