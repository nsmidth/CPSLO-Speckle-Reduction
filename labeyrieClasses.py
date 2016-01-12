# Classes used in labeyrie speckle processing

# Module Includes
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2,ifft2, fftshift
from astropy.io import fits

## Target Class: Holds data for a Reference or Binary Star
class target(): 

    fits = None				# Holds data from FITS file
    fitsFileName = None 	# Holds filename of FITS file
    psd = None				# Holds PSD of target
    psdFileName	= None		# Holds filename of PSD file
	
	# Import FITS data from filename
    #  filePath = file path of FITS file
    #  printInfo = print FITS file information?
    #  printHeaders = print FITS file headers?
    def fitsImport(self, filePath=None, printInfo = True, printHeaders = False): 
        # Use fitsFileName as filePath by default
        if filePath == None:
            filePath = self.fitsFileName				
		# Open FITS Data
        HDUList = fits.open(filePath)  
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

	# View an image from original FITS
    #  imNum = index of image to be printed/returned
    def fitsView(self, imNum=0):
        plt.figure()        
        plt.imshow(self.fits[imNum])
        plt.title('FITS Image ' + str(imNum))
        plt.show()

	# fitsExport(): Export FITS file of a np array

	# psdImport(): Import PSD data (from a FITS file)
	# psdCalc(): Calculate PSD of target
	# psdView(): View PSD


class deconvolved():
## Deconvolved Class: Holds data for devonvolved targets
# psd: Holds deconvolved PSD
# psdDeconvolve(binary, reference): Deconvolve PSDs 
# psdView(): View PSD

# psdFiltered: Holds filtered PSD
# psdFilter(LPF, HPF, Interference): Filters PSD to make PSDFiltered
# psdFilteredView(): View Filtered PSD

# acorr: Holds autocorrelation of PSD
# acorrCalc: Calculate autocorrelation from filtered PSD
# acorrView(): View autocorrelation
	pass

## Need another class to hold all this centroid finding stuff
# acorrMarked: Autocorrelation with user markings
# acorrMarkedBlank(): Erase Acorr
# acorrMark(x,y,shape): Mark acorr at location with indicated shape
# centroid[0,1][x,y]: Locations of centroids
# centroidExpected[x,y]: Expected location of secondary
# centroidEstimate(): Estimate centroids of image
# centroidCalculate(x,y,radius): Calculate centroid within area



