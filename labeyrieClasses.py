# Classes used in labeyrie speckle processing

## Target Class: Holds data for a Reference or Binary Star
class target(): 
# fits: Holds data from FITS file
# fitsFilename: Holds filename of FITS file
# fitsImport(): Import FITS data from filename
# fitsView(): View an image from original FITS

# psd: Holds PSD of target
# psdFilename: Holds filename of PSD file
# psdImport(): Import PSD data (from a FITS file)
# psdCalc(): Calculate PSD of target
# psdView(): View PSD
	pass

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



