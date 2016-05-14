# Includes
import matplotlib.pyplot as plt
import numpy as np
from classes_atmos_sim import atmospheric_simulation
from scipy.signal import argrelextrema

# Size of simulated image
nxy = 512

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
sim.phi = 0 # Set angle in degrees  

# Initializing Empty Images
sim.input_img = np.zeros((nxy,nxy)) # Input binary star object
sim.aperture_screen_s = np.zeros((nxy,nxy)) #  Aperture screen
sim.pupil_screen = np.zeros((nxy,nxy)) # Total Pupil Screen
sim.atmosphere_screen = np.zeros((nxy,nxy)) # Atmosphere screen
sim.psf = np.zeros((nxy,nxy)) # PSF of aperture
sim.binary_img = np.zeros((nxy,nxy)) # Simulated Binary Image

## Run Simulation
# Create ideal binary star image
sim.create_input_image()
# Create telescope aperture image
sim.create_aperture()
# Create atmospheric distortion phase screen
# We aren't simulating atmospheric distortion, so we just set this to a constant
sim.atmosphere_screen.fill(1)
# Calculate PSF of aperture and atmospheric phase screen
sim.get_psf()
# Calculate simulated binary star image from PSF
sim.get_binary()

#Plots
colormap = "gray"
plt.figure(figsize = (14,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(sim.aperture_screen_s, cmap=colormap, extent = (-sim.X_aperture_s_meff/2,sim.X_aperture_s_meff/2,-sim.X_aperture_s_meff/2,sim.X_aperture_s_meff/2))
plt.xlabel("[m]")
plt.ylabel("[m]")
plt.title("Aperture Screen")

plt.subplot(1,2,2)
plt.imshow(np.log10(sim.psf), cmap=colormap)
plt.xlabel("[pixels]")
plt.ylabel("[pixels]")
plt.title("Aperture PSF")

plt.figure(figsize = (14,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(sim.input_img, cmap=colormap, extent = (-nxy*sim.platescale/2,nxy*sim.platescale/2,-nxy*sim.platescale/2,nxy*sim.platescale/2))
plt.xlabel("[arcsec]")
plt.ylabel("[arcsec]")
plt.title("Input Image")

plt.subplot(1,2,2)
plt.imshow(sim.binary_img, cmap=colormap)
plt.xlabel("[pixels]")
plt.ylabel("[pixels]")
plt.title("Sensor Image")

# Show center row of aperture PSF and sensor's image
plt.figure(figsize = (20,10), dpi = 100)
dx = 100 # Number of points to be plotted around the center of image
x = np.arange(sim.center-50,sim.center+50) # Points to be plotted

# Aperture PSF Plot
plt.subplot(1,2,1)
plt.plot(x,sim.psf[sim.center,x])
plt.title("Aperture PSF (Center Row)")

# Output Image Plot
plt.subplot(1,2,2)
plt.plot(x,sim.binary_img[sim.center,x])
plt.title("Output Image (Center Row)")

# Show figures
plt.show()

# List distances from center of nulls in PSF
print("Aperture PSF Null Radii")
minima_i = argrelextrema(sim.psf[sim.center], np.less)
minima_i = np.subtract(minima_i,sim.center)[0]
print_num = 5
minima_found = 0
for value in minima_i:
    if (value > 0) and (minima_found < print_num):
        print(value)
        minima_found = minima_found + 1
