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
n_exposures = 100 # Number of simulated exposures

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
psf_avg = np.zeros((nxy,nxy))
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

  # Integrate PSFs
  psf_avg += sim.psf

  # Calculate PSD of binary star image
  binary_psd = np.power(np.abs(fftshift(fft2(sim.binary_img))),2)
  # Integrate PSDs
  binary_psd_avg += binary_psd

  # Calculate PSD of PSF
  reference_psd = np.power(np.abs(fftshift(fft2(sim.psf))),2)
  # Integrate PSDs
  reference_psd_avg += reference_psd
    
# Calculate average of PSF/PSDs
psf_avg /= n_exposures
binary_psd_avg /= n_exposures
reference_psd_avg /= n_exposures

# Calculate PSDs
input_img_psd = np.power(np.abs(fftshift(fft2(sim.input_img))),2)

# Calculate Acorrs
input_img_acorr = np.abs(fftshift(ifft2(input_img_psd)))
reference_acorr_avg = np.abs(fftshift(ifft2(reference_psd_avg)))
binary_acorr_avg = np.abs(fftshift(ifft2(binary_psd_avg)))

# Deconvolve reference from binary
deconvolved_psd_avg = np.divide(binary_psd_avg, reference_psd_avg)
deconvolved_acorr_avg = np.abs(fftshift(ifft2(deconvolved_psd_avg)))

colormap = "jet"

plt.figure(figsize = (18,8), dpi = 150)
plt.subplot(1,2,1)
plt.imshow(sim.psf, cmap=colormap)
plt.title("Short Exposure PSF")
plt.subplot(1,2,2)
plt.imshow(psf_avg, cmap=colormap)
plt.title("Long Exposure PSF")

plt.figure(figsize = (16,16), dpi = 50)
plt.subplot(3,3,1)
plt.imshow(sim.emphasized_image(img=sim.input_img,circle_radius=4), cmap=colormap)
plt.title("Binary Star Object")
plt.subplot(3,3,2)
plt.imshow(sim.psf, cmap=colormap)
plt.title("Atmospheric/Aperture PSF")
plt.subplot(3,3,3)
plt.imshow(sim.binary_img, cmap=colormap)
plt.title("Binary Star Image")
plt.subplot(3,3,4)
plt.imshow(np.log10(input_img_psd), cmap=colormap)
plt.title("Binary Star Object PSD")
plt.subplot(3,3,5)
plt.imshow(np.log10(reference_psd_avg), cmap=colormap)
plt.title("Avg PSF PSD")
plt.subplot(3,3,6)
plt.imshow(np.log10(binary_psd_avg), cmap=colormap)
plt.title("Avg Binary Star Image PSD")
plt.subplot(3,3,7)
plt.imshow(sim.emphasized_image(img=input_img_acorr,circle_radius=4), cmap=colormap)
plt.title("Binary Star Object Autocorrelation")
plt.subplot(3,3,8)
plt.imshow(reference_acorr_avg, cmap=colormap)
plt.title("Avg PSF Autocorrelation")
plt.subplot(3,3,9)
plt.imshow(binary_acorr_avg, cmap=colormap)
plt.title("Avg Binary Star Image Autocorrelation")

fig = plt.figure(figsize = (10,10), dpi = 100)
ax = fig.gca(projection='3d')
[xx,yy] = np.meshgrid(np.arange(nxy), np.arange(nxy))
ax.plot_surface(X=xx,
                Y=yy,
                Z=np.log10(reference_psd_avg),
                cmap="jet",
                linewidth=0, 
                antialiased=True)
plt.title("Reference Star PSD")

plt.figure(figsize = (14,14), dpi = 100)
plt.subplot(2,2,1)
plt.imshow(np.log10(binary_psd_avg), cmap=colormap)
plt.title("Avg Binary Star Image PSD")
plt.subplot(2,2,2)
plt.imshow( np.log10(deconvolved_psd_avg), cmap=colormap)
plt.title("Avg Deconvolved PSD")
plt.subplot(2,2,3)
plt.imshow(binary_acorr_avg, cmap=colormap)
plt.title("Avg Binary Star Image Autocorrelation")
plt.subplot(2,2,4)
plt.imshow(sim.emphasized_image(img=deconvolved_acorr_avg,circle_radius=4), cmap=colormap)
plt.title("Avg Deconvolved Autocorrelation")

plt.show()