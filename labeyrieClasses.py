# Classes used in labeyrie speckle processing

# Module Includes
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from astropy.io import fits
import sys, os

# Class to hold an array of data
# Built in methods for importing/exporting/viewing the data
class fitsData():
    # Init
    def __init__(self):
        self.data = None # Holds data
        self.fileName = None # Holds filename for import/export

    # Read data from FITS file
    # Enter the number of dimensions of the FITS file to check for it on opening
    def read(self, numDimensions=2, printInfo=True, printHeaders=False): # Imports FITS file
        # Check if input file is .fits
        if (os.path.splitext(self.fileName)[1] != ".fits"):
            # Exit program if not FITS
            sys.exit(("ERROR: " + self.fileName + " is not .fits"))

        # Open FITS Data
        HDUList = fits.open(self.fileName)
        # Print FITS File Info & Headers
        if (printInfo == True):
            print(HDUList.info())
        if (printHeaders == True):
            print("Headers:")
            print(repr(HDUList[0].header))

        # Check that input psd FITS is appropriate dimension
        if (len(np.shape(HDUList[0].data)) != numDimensions):
            sys.exit(("ERROR: " + self.fileName + " dimensions != " + str(numDimensions)))

        # Save data in FITS cube to class's data variable, then close FITS file
        self.data = HDUList[0].data
        HDUList.close()

    # Write data to FITS file
    def write(self): # Write FITS file
        # Create PrimaryHDU object with data
        hdu = fits.PrimaryHDU(self.data)
        # Create HDUList object w/ PrimaryHDU
        hdulist = fits.HDUList([hdu])
        # Write to new file
        hdulist.writeto(self.fileName)

    # View data
    # Option to view log of data
    # Option to select which image in 3D cube to view
    def view(self, log=False, title=None, imNum=0):
        plt.figure()
        # Check for dimensionality of data
        # If 2D, show the image
        if len(np.shape(self.data)) == 2:
            if (log == False):
                plt.imshow(self.data)
            if (log == True):
                plt.imshow(np.log10(self.data))
        # If 3D, show the image of selected index
        elif len(np.shape(self.data)) == 3:
            if (log == False):
                plt.imshow(self.data[imNum])
            if (log == True):
                plt.imshow(np.log10(self.data[imNum]))
        # If other, must be an error
        else:
            sys.exit("Can only view 2D or 3D data")
        plt.title(title)
        plt.show()

# Target Class: Holds data for a Reference or Binary Star
class target():
    # Init
    def __init__(self):
        self.fits = fitsData()		# Holds raw astronomical data
        self.psd = fitsData()    	# Holds calculated PSD of target

    def psdCalc(self):
    # Calculate PSD of FITS data
        # Checking if FITS data is an array of images
        if (len(self.fits.shape) == 3):
            # Generate empty array the size of an image to be used to accumulate
            #  PSD values before averaging.
            psdSum = np.zeros(self.fits.shape[1:3])

            imgNum = np.shape(self.fits)[0] # Number of images
            imgIncrement = 50 # How often to display a status message

            # Looping through all images in cube
            for index,img in enumerate(self.fits):

                # Print current file being processed
                if (((index+1) % imgIncrement) == 0):
                    print("Processed Image #: ",(index+1),"/",imgNum)

                # FFT function requires little-endian data, so casting it
                img = img.astype(float)

                # Calculate 2D power spectrum
                # This gives us only real values
                psdImg = np.abs(fft2(img))**2

                # Accumulate current PSD value
                psdSum = np.add(psdSum,psdImg)

            # Divide by # of images to calculate average
            psdAvg = np.divide(psdSum,imgNum)

            # Normalizing FFT
            psdAvg = np.divide(psdAvg, (psdAvg.size)**2)

        #Otherwise if FITS data is only one image
        elif (len(self.fits.shape) == 2):
            # FFT function requires little-endian data, so casting it
            img = self.fits.astype(float)

            # Calculate 2D power spectrum
            # This gives us only real values
            psdImg = np.abs(fft2(img))**2

            # Normalizing FFT
            psdAvg = np.divide(psdImg, (psdImg.size)**2)

        self.psd = fftshift(psdAvg)

## Deconvolved Class: Holds data for devonvolved targets
class deconvolved():
    # Init
    def __init__(self):
        self.psd = fitsData()          # Holds deconvolved PSD
        self.psdFiltered = fitsData()  # Holds filtered PSD
        self.acorr = fitsData()        # Holds autocorrelation of PSD

    # Deconvolve PSDs
    def psdDeconvolve(self, psdBinary, psdReference, constant):
        # Divide double PSD by modified reference single PSD to deconvolve
        self.psd.data = np.divide(psdBinary, (psdReference+constant))

    # psdFilter(LPF, HPF, Interference): Filters PSD to make PSDFiltered
    # LPF and HPF are gaussian shape
    # Interference filter only works for images with even number of pixels to a side
    def psdFilter(self, lpfRadius=None, hpfRadius=None, interference=False):

        imgSize = np.shape(self.psd.data)[0] # Calculate dimension of image
        imgCenter = imgSize/2         # Center index of image

        # Start filtered output as unfiltered psd
        self.psdFiltered.data = np.array(self.psd.data)

        # Perform LPF filtering if user selected a radius
        if (lpfRadius != None):

            # Save present input values
            psdTemp = np.array(self.psdFiltered.data)

            # Loop through all Y
            for y in np.arange(0,imgCenter+1):
                # Loop through all X
                for x in np.arange(0, imgCenter+1):
                    # Calculate distance of point to center of image
                    radius = np.linalg.norm((y-imgCenter,x-imgCenter))
                    # Calculate filter's value at this radius
                    filterH = np.exp(-(np.power(radius,2)/(2*np.power(lpfRadius,2))))
                    # Multiply filter's value by this point and the mirrored locations in image
                    self.psdFiltered.data[y,x] = np.multiply(psdTemp[y,x],filterH)
                    self.psdFiltered.data[y,imgSize-1-x] = np.multiply(psdTemp[y,imgSize-1-x],filterH)
                    self.psdFiltered.data[imgSize-1-y,x] = np.multiply(psdTemp[imgSize-1-y,x],filterH)
                    self.psdFiltered.data[imgSize-1-y,imgSize-1-x] = np.multiply(psdTemp[imgSize-1-y,imgSize-1-x],filterH)

        # Perform HPF filtering if user selected a radius
        if (hpfRadius != None):
            print("HPF Not Implemented Yet")
            pass

        # Perform interference filtering if user enabled it
        if (interference == True):
            # Copy image to be filtered
            psdTemp = np.array(self.psdFiltered.data)

            # Perform filtering on vertical interference
            for y in np.arange(imgSize):
                avg = np.average((psdTemp[y,imgCenter+1],psdTemp[y,imgCenter-1]))
                self.psdFiltered.data[y,imgCenter]=avg

            # Perform filtering on horizontal interference
            for x in np.arange(imgSize):
                avg = np.average((psdTemp[imgCenter+1,x],psdTemp[imgCenter-1,x]))
                self.psdFiltered.data[imgCenter,x]=avg

    # Calculate autocorrelation from filtered PSD
    def acorrCalc(self):
        # Because we normalized after FFT by multiplying by 1/N^2, and ifft
        #  function does this as well, we need to multiply by N^2 before ifft
        #  to prevent performing normalization twice
        self.acorr.data = self.psdFiltered.data*(self.psdFiltered.data.size)

        # Do iFFT on PSD's, bringing back to spatial domain
        # This should give us the autocorrelations of original images
        self.acorr.data = ifft2(self.acorr.data)

        # Taking iFFT of PSD (all real values) results in complex valued output
        #  Must view the magnitude of the output
        #  Doing FFTshift to move eyes of autocorrelation near center
        # Taking iFFT of PSD (all real values) results in complex valued output
        #  Must view the magnitude of the output
        self.acorr.data = np.abs(fftshift(self.acorr.data))



## Need another class to hold all this centroid finding stuff
# acorrMarked: Autocorrelation with user markings
# acorrMarkedBlank(): Erase Acorr
# acorrMark(x,y,shape): Mark acorr at location with indicated shape
# centroid[0,1][x,y]: Locations of centroids
# centroidExpected[x,y]: Expected location of secondary
# centroidEstimate(): Estimate centroids of image
# centroidCalculate(x,y,radius): Calculate centroid within area

