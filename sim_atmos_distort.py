# Includes
import matplotlib.pyplot as plt
import numpy as np
from classes_atmos_sim import atmospheric_simulation

nxy = 512

sim = atmospheric_simulation()

sim.create_input_image()
sim.create_aperture()
sim.create_atmosphere_screen()
sim.get_psf()
sim.get_binary()
img_noisy = sim.add_noise(sim.psf,photons=1E5,gaussian_var=1E-3)

## Adding photon shot noise
photons = (1000,100000)
photons_n = len(photons)

img_shot_noise = np.zeros((photons_n,nxy,nxy))

for i in np.arange(photons_n):
  # Use image as lambda to get random values based on poisson distribution at each point
  img_shot_noise[i] = sim.add_noise(sim.binary_img,photons=photons[i],gaussian_var=0)

##Plots
colormap = "jet"

plt.figure(figsize = (18,8))
plt.subplot(1,2,1)
plt.imshow(sim.phase_screen, cmap=colormap)
plt.colorbar()
plt.title("Atmosphere Phase Shift (unwrapped)")
plt.subplot(1,2,2)
plt.imshow(sim.aperture_screen_s, cmap=colormap)
plt.title("Telescope Aperture")

plt.figure(figsize = (8,8))
plt.imshow(sim.psf, cmap=colormap)
plt.title("Total System PSF")

plt.figure(figsize = (14,6))
plt.subplot(1,2,1)
plt.imshow(sim.emphasized_image(sim.input_img), cmap=colormap, extent = (-6.513,6.513,-6.513,6.513))
plt.xlabel("[arcsec]")
plt.ylabel("[arcsec]")
plt.title("Input Image")
plt.subplot(1,2,2)
plt.imshow(sim.binary_img, cmap=colormap)
plt.xlabel("[pixels]")
plt.ylabel("[pixels]")
plt.title("Binary Image")

plt.figure(figsize = (10+10*photons_n,10), dpi = 50)
plt.subplot(1,1+photons_n,1)
plt.imshow(sim.binary_img, cmap=colormap)
plt.title("Sensor Image (No Shot Noise)")
for i in np.arange(2):
  plt.subplot(1,1+photons_n,2+i)
  plt.imshow(img_shot_noise[i])
  plt.title(("Sensor Image (" + str(photons[i]) + " Photons)"))

  
plt.show()



