# Module for performing full Labeyrie deconvolution of FITS files

# Included modules
import sys
from labeyrieClasses import target,deconvolved

# Instantiate objects
binary = target()
reference = target()
dcnvlv = deconvolved()

# Prompt user if they'd like to preprocess:
if ( input("Preprocessing raw FITS files? [y/n]   ").lower() == 'y' ):
	# Prompt for binary FITS file location
	#binary.fitsFileName = 
	# Import binary FITS
	#binary.fitsImport()
	# Preprocess
	#binary.PSDCalc()

	# Prompt for reference FITS file location
	#reference.fitsFileName = 
	# Import reference FITS
	#reference.fitsImport()
	# Preprocess
	#reference.PSDCalc()
	
	pass
else:
	# Prompt user for binary PSD file location
	#binary.psdFileName = 
	# Import binary star PSD FITS
	#binary.psdImport()
	
	# Prompt user for reference PSD file location
	#reference.psdFileName = 
	# Import reference star PSD FITS
	#reference.psdImport()
	pass
	
# Deconvolve reference and binary stars
#dcnvlv.psdDeconvolve(binary.PSD,reference.PSD)

# Perform filtering on output star object
#dcnvlv.psdFilter(lpf = 15, interference = TRUE)

# Get autocorrelogram
#dcnvlv.acorrCalc()

# Return autocorrelogram as PNG or FITS
