# Classes used in labeyrie speckle processing

# Module Includes
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from astropy.io import fits
import sys, os


# Target Class: Holds data for a Reference or Binary Star
class target():

    fits = None				# Holds data from FITS file
    fitsFileName = None 	# Holds filename of FITS file
    psd = None				# Holds PSD of target
    psdFileName	= None		# Holds filename of PSD file

    # Import FITS data from filename
    #  filePath = file path of FITS file
    #  printInfo = print FITS file information?
    #  printHeaders = print FITS file headers?
    def fitsImport(self, printInfo = True, printHeaders = False):

        # Check if input file is .fits
        if (os.path.splitext(self.fitsFileName)[1] != ".fits"):
            # Exit program if not FITS
            sys.exit(("ERROR: " + self.fitsFileName + " is not .fits"))

        # Open FITS Data
        HDUList = fits.open(self.fitsFileName)
        # Print FITS File Info & Headers
        if (printInfo == True):
            print(HDUList.info())
        if (printHeaders == True):
            print("Headers:")
            print(repr(HDUList[0].header))
        # Save data in FITS cube to local variable, then close FITS file
        fitsData = HDUList[0].data
        HDUList.close()
        # Save fits data
        self.fits = fitsData

    # View an image from original FITS
    #  imNum = index of image to be printed/returned
    def fitsView(self, imNum=0):
        plt.figure()
        plt.imshow(self.fits[imNum])
        plt.title('FITS Image ' + str(imNum))
        plt.show()

    ## fitsExport(): Export FITS file of fits data
    def fitsExport(self):

        # Create PrimaryHDU object with data
        hdu = fits.PrimaryHDU(self.fits)
        # Create HDUList object w/ PrimaryHDU
        hdulist = fits.HDUList([hdu])
        # Write to new file
        hdulist.writeto(self.fitsFileName)

    ## psdImport(): Import PSD data (from a FITS file)
    def psdImport(self, printInfo = True, printHeaders = False):

        # Check if input file is .fits
        if (os.path.splitext(self.psdFileName)[1] != ".fits"):
            # Exit program if not FITS
            sys.exit(("ERROR: " + self.psdFileName + " is not .fits"))

        # Open FITS Data
        HDUList = fits.open(self.psdFileName)
        # Print FITS File Info & Headers
        if (printInfo == True):
            print(HDUList.info())
        if (printHeaders == True):
            print("Headers:")
            print(repr(HDUList[0].header))

        # Check that input psd FITS is 2D
        if (len(np.shape(HDUList[0].data)) > 2):
            sys.exit(("ERROR: " + self.psdFileName + " contains more than one image"))

        # Save data in FITS cube to local variable, then close FITS file
        psdData = HDUList[0].data
        HDUList.close()
        # Save fits data
        self.psd = psdData

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

    # psdView(): View PSD
    def psdView(self):
        plt.figure()
        plt.imshow(np.log10(self.psd))
        plt.title('PSD Image ')
        plt.show()

    ## psdExport(): Export FITS file of psd
    def psdExport(self):
        # Create PrimaryHDU object with data
        hdu = fits.PrimaryHDU(self.psd)
        # Create HDUList object w/ PrimaryHDU
        hdulist = fits.HDUList([hdu])
        # Write to new file
        hdulist.writeto(self.psdFileName)

## Deconvolved Class: Holds data for devonvolved targets
class deconvolved():
    psd = None          # Holds deconvolved PSD
    psdFiltered = None  # Holds filtered PSD
    acorr = None        # Holds autocorrelation of PSD

    # Deconvolve PSDs
    def psdDeconvolve(self, psdBinary, psdReference, constant):
        # Divide double PSD by modified reference single PSD to deconvolve
        self.psd = np.divide(psdBinary, (psdReference+constant))

    # View PSD
    def psdView(self):
        plt.figure()
        plt.imshow(np.log10(self.psd))
        plt.title('PSD Image')
        plt.show()

# psdFiltered: Holds filtered PSD
# psdFilter(LPF, HPF, Interference): Filters PSD to make PSDFiltered
# psdFilteredView(): View Filtered PSD

# acorr: Holds autocorrelation of PSD
# acorrCalc: Calculate autocorrelation from filtered PSD
# acorrView(): View autocorrelation


## Need another class to hold all this centroid finding stuff
# acorrMarked: Autocorrelation with user markings
# acorrMarkedBlank(): Erase Acorr
# acorrMark(x,y,shape): Mark acorr at location with indicated shape
# centroid[0,1][x,y]: Locations of centroids
# centroidExpected[x,y]: Expected location of secondary
# centroidEstimate(): Estimate centroids of image
# centroidCalculate(x,y,radius): Calculate centroid within area



