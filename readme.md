# CPSLO-Speckle-Reduction

Creating a program for reduction of binary star FITS images for astrometric and photometric data

- dev/ : All experimental data/iPython notebooks go here. Algorithms often developed here first
- labeyrie.py : Performs Labeyrie deconvolution on the average PSD of a binary stars/reference star
- preprocess.py : Calculates average PSD of FITS files. Saves results in _PSD.fits files
- labeyrieClasses : Holds all classes/functions used in these scripts

Starting from scratch, one would reduce speckle data by using the following procedure:

1. Run preprocess.py script.
2. Select all .fits speckle data and run. This produces _PSD.fits files in the same folder as the original data
3. Run labeyrie.py script
4. Select a binary/reference star pair and run. 
5. Autocorrelation is calculated and displayed
6. User is prompted for filename for saved autocorrelation data (X out the dialog if you don't want to save)
7. User is prompted if they'd like to do a centroid estimation of secondary location
7. User must enter photometric data, then calculations are done and shown
