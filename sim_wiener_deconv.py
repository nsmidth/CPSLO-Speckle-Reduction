# Includes
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.signal import fftconvolve
from scipy.signal import argrelextrema
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from classes_atmos_sim import atmospheric_simulation

## Parameters
# Image specs
nxy = 512
center = int(nxy/2)
n_exposures = 1 # Number of simulated exposures

# Instantiate simulation object
sim = atmospheric_simulation()

## Assign all variables of simulation
# Aperture/Telescope Specifications:
sim.diameter_m = 2.133 # Mirror diameter in meters
sim.focal_length = 129.69 # Effective focal length in meters

# Camera Specs
sim.wavelength = 0.8E-6 # Wavelength of light
sim.pixel = 8E-6 # length of pixel side in m
sim.bits = 14 # Bits in camera ADC
sim.nxy = nxy # Length of sensor side in pixels
sim.center = int(nxy/2) # Center of sensor in pixels
sim.platescale = 206265*sim.pixel/(sim.focal_length) # Calculate plate scale
sim.gamma = 1.6 # Gamma correction

# Binary star specs
sim.rho = 0.5 # Set separation in arcseconds
sim.phi = 45 # Set angle in degrees  

# Atmospheric Specs
sim.alpha = 1/100 # Multiplicative constant
sim.r0 = 0.2 # Fried Parameter

# Initializing Empty Images
sim.input_img = np.zeros((nxy,nxy)) # Input binary star object
sim.aperture_screen_s = np.zeros((nxy,nxy)) #  Aperture screen
sim.phase_screen = np.zeros((nxy,nxy)) #  Phase screen values (unwrapped)
sim.atmosphere_screen = np.zeros((nxy,nxy)) #  Atmosphere screen
sim.pupil_screen = np.zeros((nxy,nxy)) # Total Pupil Screen
sim.psf = np.zeros((nxy,nxy)) # PSF of aperture/atmosphere
sim.binary_img = np.zeros((nxy,nxy)) # Simulated Binary Image

## Run Simulation
# Create image accumulation arrays
binary_psd_avg = np.zeros((nxy,nxy))
reference_psd_avg = np.zeros((nxy,nxy))

# Create ideal binary star image
sim.create_input_image()
# Create telescope aperture image
sim.create_aperture()
# Generate multiple images, integrate them
for i in np.arange(n_exposures):
  # Create atmospheric phase screen image
  sim.create_atmosphere_screen()
  # Calculate PSF of aperture and atmospheric phase screen
  sim.get_psf()
  # Calculate simulated binary star image from PSF
  sim.get_binary()
  # Add noise to binary and reference stars
  sim.psf = sim.add_noise(sim.psf,1E5,1E-4)
  sim.binary_img = sim.add_noise(sim.binary_img,1E5,1E-4)
  
  # Calculate PSD of binary star image
  binary_psd = np.power(np.abs(fftshift(fft2(sim.binary_img))),2)
  # Integrate PSDs
  binary_psd_avg += binary_psd

  # Calculate PSD of PSF
  reference_psd = np.power(np.abs(fftshift(fft2(sim.psf))),2)
  # Integrate PSDs
  reference_psd_avg += reference_psd
    
# Calculate average of PSF/PSDs
binary_psd_avg /= n_exposures
reference_psd_avg /= n_exposures

# Calculate Acorrs
reference_acorr_avg = np.abs(fftshift(ifft2(reference_psd_avg)))
binary_acorr_avg = np.abs(fftshift(ifft2(binary_psd_avg)))

# Deconvolve reference from binary
deconvolved_psd_avg = np.divide(binary_psd_avg, reference_psd_avg)

# Calculate Autocorr
deconvolved_acorr_avg = np.abs(fftshift(ifft2(deconvolved_psd_avg)))

colormap = "jet"

plt.figure(figsize = (14,14), dpi = 100)
plt.subplot(2,3,1)
plt.imshow(np.log10(reference_psd_avg), cmap=colormap)
plt.title("Reference Star Image PSD")
plt.subplot(2,3,2)
plt.imshow(np.log10(binary_psd_avg), cmap=colormap)
plt.title("Binary Star Image PSD")
plt.subplot(2,3,3)
plt.imshow( np.log10(deconvolved_psd_avg), cmap=colormap)
plt.title("Deconvolved PSD")
plt.subplot(2,3,4)
plt.imshow(reference_acorr_avg, cmap=colormap)
plt.title("Reference Star Image PSD")
plt.subplot(2,3,5)
plt.imshow(binary_acorr_avg, cmap=colormap)
plt.title("Binary Star Image Autocorrelation")
plt.subplot(2,3,6)
plt.imshow(sim.emphasized_image(deconvolved_acorr_avg), cmap=colormap)
plt.title("Deconvolved Autocorrelation")

plt.show()