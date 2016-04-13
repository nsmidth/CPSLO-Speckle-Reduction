# Classes used in labeyrie speckle processing

# Module Includes
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from astropy.io import fits
import sys, os, cv2

# Class to hold an array of data
# Built in methods for importing/exporting/viewing the data
# Used in target and deconvolved class
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
        plt.colorbar()
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
        if (len(self.fits.data.shape) == 3):
            # Generate empty array the size of an image to be used to accumulate
            #  PSD values before averaging.
            psdSum = np.zeros(self.fits.data.shape[1:3])

            imgNum = np.shape(self.fits.data)[0] # Number of images
            imgIncrement = imgNum/20 # How often to display a status message

            # Looping through all images in cube
            for index,img in enumerate(self.fits.data):

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

        self.psd.data = fftshift(psdAvg)

# Deconvolved Class: Holds data for devonvolved targets
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

            # Create centered meshgrid of image
            xx,yy = np.meshgrid(np.arange(imgSize),np.arange(imgSize))
            xx = np.subtract(xx,imgCenter)
            yy = np.subtract(yy,imgCenter)
            rr = np.power(np.power(xx,2)+np.power(yy,2),0.5)

            # Create LPF filter image
            filterH = np.exp(-(np.power(rr,2)/(2*np.power(lpfRadius,2))))

            # Filter the PSD with LPF
            self.psdFiltered.data = np.multiply(psdTemp,filterH)

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

# coords class: Holds methods for working in polar/cart cordinates
# Used in camsky class
class coords():
    def __init__(self, midpoint):
        self.theta = None
        self.rho = None
        self.x = None
        self.y = None
        self.midpoint = midpoint

    # Orbital plot has 0deg pointing down on the screen [resulting in 90deg shift]
    # Theta rotates counterclockwise from 0deg
    # Positive pixel direction is downward [resulting in Y phase reversal]
    # Polar coordinates centered on middle pixel of image

    # convert theta and rho values to image x and y coordinates
    def polar2cart(self):
        self.x = self.midpoint + (self.rho * np.cos(np.radians(self.theta-90)))
        self.y = self.midpoint - (self.rho * np.sin(np.radians(self.theta-90)))

    # convert x and y coordinates to theta and rho values
    def cart2polar(self):
        # Get cartesians centered on 0 for calculations
        x = self.x-self.midpoint
        y = self.midpoint-self.y

        self.rho = np.linalg.norm((x,y))
        self.theta = np.degrees(np.arctan2(y,x))+90
        # Add 360 deg if theta is a negative angle
        self.theta = self.theta + 360*(self.theta<0)

# camsky class: Holds methods for working between measurements in camera and sky
class camsky():
    def __init__(self, midpoint, delta, e):
        self.cam = coords(midpoint)
        self.sky = coords(midpoint)
        self.delta = delta
        self.e = e

    # Convert from polar coords in sky to polar coords in camera
    def sky2cam(self):
        self.cam.theta = self.sky.theta+self.delta
        self.cam.rho = self.sky.rho/self.e

    # Convert from polar coords in camera to polar coords in sky
    def cam2sky(self, camTheta, camRho):
        self.sky.theta = self.cam.theta-self.delta
        self.sky.rho = self.cam.rho*self.e

# Photometry class: Holds data for photometry measurements
class photometry():
    # Init
    def __init__(self):
        self.acorr = None # Raw Autocorrelation
        self.acorrMarked = None # Marked up Autocorrelation

    # Clear marked up acorr back to original acorr
    def acorrMarkedClear(self):

        # Create blank RGB image
        imgDim = np.shape(self.acorr)[0]
        self.acorrMarked = np.zeros((imgDim,imgDim,3), np.uint8)

        # Calculate scaled version of autocorr image
        imgScaled = np.divide(self.acorr,np.max(self.acorr))
        imgScaled = (np.multiply(imgScaled, 255)).astype(np.uint8)

        # Fill in scaled image into each R/G/B Channel
        for i in np.arange(3):
            self.acorrMarked[:,:,i] = imgScaled

    # Mark acorr at location with indicated shape
    def acorrMark(self,x,y,shape,color,radius=10):
        if shape == 'o':
            cv2.circle(self.acorrMarked, (x, y), radius, color, 1)
        elif shape == '+':
            cv2.line(self.acorrMarked,(x-radius,y),(x+radius,y),color,1)
            cv2.line(self.acorrMarked,(x,y-radius),(x,y+radius),color,1)
        else:
            print("Invalid shape input")

# centroidExpected[x,y]: Expected location of secondary
# centroidCalculate(x,y,radius): Calculate centroid within area

    # Estimate centroids of autocorr
    # Autocorrelation image consists primarily of the large central lobe and less tall two side lobes
    # Lobes are found by creating a thresholded image, 1 where above threshold otherwise 0
    # Threshold is started at peak of center lobe, moved down to the minimum value that gives 3 contours
    def centroidEstimate(self):

        # Start thresh at max value of autocorr (assuming that is from central peak)
        threshold = np.floor(self.acorr.max())

        # Calculate values to increment/decrement threshold by depending on height of autocorrelogram
        incrementFine   = self.acorr.max()/100

        # Assume we start with 1 contour visible
        numContours = 1

        # Step down threshold until we see the 3 contours we expect
        while (numContours != 3):
            # Decrement the threshold
            threshold = threshold - incrementFine

            # Create thresholded image
            acorrThresh = self.acorr.data > threshold # Calculate indices in image above threshold
            acorrThresh = (np.multiply(acorrThresh,255).astype(np.uint8))

            # Find contours of threshold image
            im,contours,hierarchy = cv2.findContours(acorrThresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
            numContours = len(contours)

            if threshold < 0:
                print("Error in finding 3 contours")
                sys.exit()

        # Now we've found our 3 main contours. Keep stepping down until we don't see 3 anymore
        # Can be less than 3 if they start to blend together, or more than 3 if other noise in
        #  image appears
        while (numContours == 3):
            # Decrement the threshold
            threshold = threshold - incrementFine

            # Create thresholded image
            acorrThresh = self.acorr.data > threshold # Calculate indices in image above threshold
            acorrThresh = (np.multiply(acorrThresh,255).astype(np.uint8))

            # Find contours of threshold image
            im,contours,hierarchy = cv2.findContours(acorrThresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
            numContours = len(contours)

            if threshold < 0:
                print("Error in moving threshold below 3 contours")
                sys.exit()

        # We've moved just below our desired 3 contours threshold. We moved down by incrementFine, so we
        #  need to move up by incrementFine to get back to the proper threshold
        threshold = threshold + incrementFine
        # Create thresholded image
        acorrThresh = self.acorr.data > threshold # Calculate indices in image above threshold
        acorrThresh = (np.multiply(acorrThresh,255).astype(np.uint8))
        # Find contours of threshold image
        im,contours,hierarchy = cv2.findContours(acorrThresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        numContours = len(contours)

        # Now that we have the contours, we want to find their centroids for sorting purposes
        centroid = np.zeros((3,2))
        centerDistance = np.zeros(3)
        imgCenter = np.divide(np.shape(self.acorr),2) # Center pixel indices

        for i in np.arange(3):
            # Calculating Image Moment/Raw Moment of contour
            M = cv2.moments(contours[i])

            # Calculate centroid X and Y components
            cx = int(M['m10']/M['m00'])
            cy = int(M['m01']/M['m00'])

            # Save centroid
            centroid[i] = (cx, cy)

            # Calculate distance from centroid to image center
            centerDistance[i] = np.linalg.norm(centroid[i] - imgCenter)

        # Get index of central contour as contour with min distance to center
        iCenter = np.where(centerDistance == centerDistance.min())
        iCenter = iCenter[0]

        # Finding indices of side lobes
        iSide = [0,1,2] # Possible indices of side lobes
        iSide.remove(iCenter) # Remove the center lobe index

        # Make empty mask images
        maskSide = np.zeros((2,np.shape(self.acorr)[0],np.shape(self.acorr)[1]))
        lobeSide = np.zeros((2,np.shape(self.acorr)[0],np.shape(self.acorr)[1]))

        # Create masks and calculate side lobes
        for i in np.arange(2):
            cv2.drawContours(maskSide[i],contours,iSide[i],(1,1,1),-1)
            lobeSide[i] = (np.multiply(self.acorr,maskSide[i]))

        # Calculating intensity weighted centroid of side lobes, using entire contour
        M00 = np.zeros((2))
        M10 = np.zeros((2))
        M01 = np.zeros((2))
        centroid = np.zeros((2,2))

        # Calculate for each side lobe
        for lobe in np.arange(2):

            # Calculate average
            M00[lobe] = np.sum(lobeSide[lobe])

            # Calculate X component
            # Loop through each column
            for column in np.arange(np.shape(lobeSide[lobe])[1]):
                M10[lobe] += np.multiply(column,np.sum(lobeSide[lobe,:,column]))

            # Calculate Y component
            # Loop through each row
            for row in np.arange(np.shape(lobeSide[lobe])[0]):
                M01[lobe] += np.multiply(row,np.sum(lobeSide[lobe,row,:]))

            # Calculate X and Y Centroid Components of each lobe
            centroid[lobe,0] = M10[lobe]/M00[lobe] # X Component
            centroid[lobe,1] = M01[lobe]/M00[lobe] # Y Component

        # Return centroid locations
        return centroid
