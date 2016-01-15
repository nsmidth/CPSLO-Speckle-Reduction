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


    # psdFilter(LPF, HPF, Interference): Filters PSD to make PSDFiltered
    # LPF and HPF are gaussian shape
    # Interference filter only works for images with even number of pixels to a side
    def psdFilter(self, lpfRadius=None, hpfRadius=None, interference=False):

        imgSize = np.shape(self.psd)[0] # Calculate dimension of image
        imgCenter = imgSize/2         # Center index of image

        # Start filtered output as unfiltered psd
        self.psdFiltered = np.array(self.psd)

        # Perform LPF filtering if user selected a radius
        if (lpfRadius != None):

            # Save present input values
            psdTemp = np.array(self.psdFiltered)

            # Loop through all Y
            for y in np.arange(0,imgCenter+1):
                # Loop through all X
                for x in np.arange(0, imgCenter+1):
                    # Calculate distance of point to center of image
                    radius = np.linalg.norm((y-imgCenter,x-imgCenter))
                    # Calculate filter's value at this radius
                    filterH = np.exp(-(np.power(radius,2)/(2*np.power(lpfRadius,2))))
                    # Multiply filter's value by this point and the mirrored locations in image
                    self.psdFiltered[y,x] = np.multiply(psdTemp[y,x],filterH)
                    self.psdFiltered[y,imgSize-1-x] = np.multiply(psdTemp[y,imgSize-1-x],filterH)
                    self.psdFiltered[imgSize-1-y,x] = np.multiply(psdTemp[imgSize-1-y,x],filterH)
                    self.psdFiltered[imgSize-1-y,imgSize-1-x] = np.multiply(psdTemp[imgSize-1-y,imgSize-1-x],filterH)

        # Perform HPF filtering if user selected a radius
        if (hpfRadius != None):
            pass

        # Perform interference filtering if user enabled it
        if (interference == True):
            # Copy image to be filtered
            psdTemp = np.array(self.psdFiltered)

            # Perform filtering on vertical interference
            for y in np.arange(imgSize):
                avg = np.average((psdTemp[y,imgCenter+1],psdTemp[y,imgCenter-1]))
                self.psdFiltered[y,imgCenter]=avg

            # Perform filtering on horizontal interference
            for x in np.arange(imgSize):
                avg = np.average((psdTemp[imgCenter+1,x],psdTemp[imgCenter-1,x]))
                self.psdFiltered[imgCenter,x]=avg

            self.psdFiltered=psdTemp


    # View Filtered PSD
    def psdFilteredView(self):
        plt.figure()
        plt.imshow(np.log10(self.psdFiltered))
        plt.title('Filtered PSD Image')
        plt.show()


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



