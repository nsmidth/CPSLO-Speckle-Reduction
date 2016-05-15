# CPSLO-Speckle-Reduction

Creating a program for reduction of binary star FITS images for astrometric and photometric data

- dev/ : All experimental data/iPython notebooks go here. Algorithms often developed here first

- bench_psd_calc.py : Benchmarking script for different methods of PSD calculation
- cython_psd.pyx : Cython PSD calculation module, must be compiled to be run
- cython_psd_make.py : Used to compile a cython_psd.so shared library
- cython_psd.so : Compiled shared library
- fftw_psd.c : C function to calculate PSD
- fftw_psd_makefile : makefile to compile fftw_psd.c to shared library, to be wrapped in Python
- fftw_psd.so : Compiled shared library

- classes_labeyrie.py : Holds all classes/functions used in FITS import/export and processing speckle data
- classes_atmos_sim.py : Holds all classes/function used in atmospheric/aperture distortion simulations
- classes_astrometry.py : Experimental modules used for finding peaks in Autocorrelation images

- labeyrie_deconv.py : Performs Labeyrie deconvolution on the average PSD of a binary stars/reference star
- labeyrie_preprocess.py : Calculates average PSD of FITS files. Saves results in _PSD.fits files

- sim_atmos_distort.py : View simulated images of binary stars through atmospheric distortion
- sim_diffraction_limit.py : Simulate observations through finite aperture with no atmospheric distortion
- sim_kpno_deconv.py : Perform Labeyrie deconvolution on real FITS data
- sim_labeyrie_deconv.py : Perform Labeyrie deconvolution on noiseless simulated data
- sim_wiener_deconv.py : Perform Labeyrie deconvolution on noisy simulated data with a Wiener filter
- sim_prime_focus.py : Create ray trace diagram for the KPNO telescope

Starting from scratch, one would reduce speckle data by using the following procedure:

1. Run labeyrie_preprocess.py script.
2. Select all .fits speckle data and run. This produces _PSD.fits files in the same folder as the original data
3. Run labeyrie_deconv.py script
4. Select a binary/reference star pair and run. 
5. Autocorrelation is calculated and displayed
6. User is prompted for filename for saved autocorrelation data (X out the dialog if you don't want to save)
7. User is prompted if they'd like to do a centroid estimation of secondary location
7. User must enter photometric data, then calculations are done and shown
