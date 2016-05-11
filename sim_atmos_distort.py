# Includes
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.signal import fftconvolve
from scipy.signal import argrelextrema

## Parameters
# Image specs
nxy = 512
center = int(nxy/2)

# Aperture/Telescope Specifications:
diameter_m = 2.133 # Mirror diameter in meters
focal_length = 129.69 # Effective focal length in meters
wavelength = 0.8E-6 # Wavelength of light
pixel = 8E-6 # Dimension of a pixel
diameter_m = 2.133 # Telescope diameter in m

## Creating Binary Star Input Image
# For KP 2.1m telescope, 
# input image is in units of 25.44 milliarcseconds
platescale = 0.02544 # Plate scale in arcsec/pixel
rho = 1.2 # Set separation in arcseconds
phi = 25 # Set angle in degrees
# Calculate coordinates of stars
x = int( rho/(2*platescale) * np.cos(np.deg2rad(phi)) )
y = int( rho/(2*platescale) * np.sin(np.deg2rad(phi)) )
x1 = center + x
y1 = center + y
x2 = center - x
y2 = center - y    
# Empty input image
input_img = np.zeros((nxy,nxy)) 
# Place stars on image
input_img[y1,x1] = 1 
input_img[y2,x2] = 1
# Scale image power to 1
input_img_power = np.sum(np.power(input_img,2))
input_img = np.divide(input_img,np.sqrt(input_img_power))

## Telescope aperture creation:
# Total spatial sample range
X_aperture_s = 1/pixel 
# dx value for sampled aperture image
dx_aperture_s = X_aperture_s/nxy 
# Coordinates of sampled image
x_aperture_s = np.arange(0,X_aperture_s,dx_aperture_s) 
x_aperture_s = x_aperture_s - X_aperture_s/2
# Meshgrid of sampled coordinates
xx_s,yy_s = np.meshgrid(x_aperture_s,x_aperture_s)
# Scaled aperture diameter to effectively 
#  resample aperture image
diameter_s = diameter_m/(focal_length*wavelength)
# Draw new circle at correct dimensions
# Calculate grid of distances from circle center
circle_s = (xx_s) ** 2 + (yy_s) ** 2 
# Draw boolean circle
circle_s = circle_s < (diameter_s/2)**2 
# Convert boolean circle to int
circle_s= circle_s.astype(np.int64)
# Save aperture image in units of meters
aperture_screen_s = circle_s
# Scale aperture image power to 1
aperture_screen_power = np.sum(np.power(aperture_screen_s,2))
aperture_screen_s = np.divide(aperture_screen_s,np.sqrt(aperture_screen_power))
# Calculate effective size of sampled aperture image in meters
X_aperture_s_meff = focal_length*wavelength/pixel

## Phase screen creation:
# Total array sample size
d_aperture = X_aperture_s_meff
# Fried parameter [m]
r0 = 0.2
# Spatial sample resolution
dxy = d_aperture/nxy
# Spatial frequency resolution
df = 1/(d_aperture) 
# Image sample indices array
x = np.multiply( np.subtract(np.arange(nxy),int(nxy/2)), dxy )
# Spatial Frequency indices array
xf = np.multiply( np.subtract(np.arange(nxy),int(nxy/2)), df )
# Meshgrid of spatial frequency domain
[xx,yy]=np.meshgrid(xf,xf)
# Radius from center meshgrid
rr = (np.sqrt(np.power(xx,2)+np.power(yy,2)))
# Calculate Kolmogorov spectral density
alpha = 1/100
phase_PSD = np.power(rr,-11/3)
phase_PSD = np.multiply(alpha*0.023/(r0**(5/3)),phase_PSD)
# Set DC component to 0 (previous calc attempts to set to 1/0)
phase_PSD[int(nxy/2),int(nxy/2)] = 0 
# Generate random phase component
random_img = np.random.normal(loc=0,scale=1,size=(nxy,nxy))
phase_phase = np.multiply(np.pi,random_img)
# Construct phase screen spectrum
phase_screen_f = np.multiply(np.sqrt(phase_PSD),np.exp(1j*phase_phase))
# Calculate phase screen
phase_screen = np.real(ifft2(fftshift(phase_screen_f)*nxy*nxy))
# Create complex atmospheric screen
atmosphere_screen = np.exp(np.multiply(1j,phase_screen))
# Generate total screen, combining atmosphere and aperture
pupil_screen = np.multiply(atmosphere_screen,aperture_screen_s)

## Calculate system's total response 
# Calculate total PSF of system
psf = fftshift(fft2(pupil_screen))
psf = np.power(np.abs(psf),2)
# Normalize PSF
psf_power = np.sum(np.power(psf,2))
psf = np.divide(psf,psf_power)
# Gamma correct PSF
gamma = 1.6
psf = np.power(psf,1/gamma)
# Convolve PSF with input image using FFT
sensor_img = fftconvolve(input_img,psf)
# Save the center 512x512 image
sensor_img = sensor_img[center:center+nxy,center:center+nxy] 

## Adding photon shot noise
# Normalize image to average of 1 electron received over whole image
img_normalized = sensor_img/(np.sum(sensor_img))

# Average number of photons incident on image
photons = (1000,100000)
photons_n = len(photons)
img_shot_noise = np.zeros((photons_n,nxy,nxy))
psd_shot_noise = np.zeros((photons_n,nxy,nxy))
acorr_shot_noise = np.zeros((photons_n,nxy,nxy))

for i in np.arange(photons_n):
  # Scale image up to the desired number of photons received
  img_scaled = img_normalized*photons[i]
  # Use image as lambda to get random values based on poisson distribution at each point
  img_shot_noise[i] = np.random.poisson(lam=np.abs(img_scaled), size=None)

##Plots
colormap = "jet"

plt.figure(figsize = (18,8))
plt.subplot(1,2,1)
plt.imshow(phase_screen, cmap=colormap)
plt.colorbar()
plt.title("Atmosphere Phase Shift (unwrapped)")
plt.subplot(1,2,2)
plt.imshow(aperture_screen_s, cmap=colormap)
plt.title("Telescope Aperture")

plt.figure(figsize = (8,8))
plt.imshow(psf, cmap=colormap)
plt.title("Total System PSF")

plt.figure(figsize = (14,6))
plt.subplot(1,2,1)
plt.imshow(input_img, cmap=colormap, extent = (-6.513,6.513,-6.513,6.513))
plt.xlabel("[arcsec]")
plt.ylabel("[arcsec]")
plt.title("Input Image")
plt.subplot(1,2,2)
plt.imshow(sensor_img, cmap=colormap)
plt.xlabel("[pixels]")
plt.ylabel("[pixels]")
plt.title("Sensor Image")

plt.figure(figsize = (10+10*photons_n,10), dpi = 50)
plt.subplot(1,1+photons_n,1)
plt.imshow(sensor_img, cmap=colormap)
plt.title("Sensor Image (No Shot Noise)")
for i in np.arange(2):
  plt.subplot(1,1+photons_n,2+i)
  plt.imshow(img_shot_noise[i])
  plt.title(("Sensor Image (" + str(photons[i]) + " Photons)"))

  
plt.show()



