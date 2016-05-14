# Includes
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.signal import fftconvolve
from scipy.signal import argrelextrema

# Image specs
nxy = 512
center = int(nxy/2)

## Calculate scaling factors for sampled aperture image
diameter_m = 2.133 # Mirror diameter in [m]
pixel = 8E-6 # Dimension of a pixel
focal_length = 129.69 # Effective focal length in meters
wavelength = 0.8E-6 # Wavelength of light
platescale = 206265 * pixel / focal_length # Calculate system's plate scale

## Creating Binary Star Input Image
rho = 0.5 # Set separation in arcseconds
phi = 0 # Set angle in degrees
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

# Total spatial sample range
X_aperture_s = 1/pixel 
# dx value for sampled aperture image
dx_aperture_s = X_aperture_s/nxy 
# Coordinates of sampled image
x_aperture_s = np.arange(0,X_aperture_s,dx_aperture_s) - X_aperture_s/2
# Meshgrid of sampled coordinates
xx_s,yy_s = np.meshgrid(x_aperture_s,x_aperture_s)
# Scaled aperture diameter to effectively resample aperture image
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

## Calculate Output Images
# Calculate normalized complex diffraction pattern
aperture_psf = fftshift(fft2(aperture_screen_s))
# Calculate PSF of complex diffraction pattern
aperture_psf = np.power(np.abs(aperture_psf),2)
# Scale PSF power to 1
aperture_psf_power = np.sum(np.power(aperture_psf,2))
aperture_psf = np.divide(aperture_psf,np.sqrt(aperture_psf_power))
# Convolve PSF with input image
sensor_img = fftconvolve(input_img,aperture_psf)
# Save the center nxy by nxy image
sensor_img = sensor_img[center:center+nxy,center:center+nxy] 

#Plots
colormap = "gray"
plt.figure(figsize = (14,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(aperture_screen_s, cmap=colormap, extent = (-X_aperture_s_meff/2,X_aperture_s_meff/2,-X_aperture_s_meff/2,X_aperture_s_meff/2))
plt.xlabel("[m]")
plt.ylabel("[m]")
plt.title("Aperture Screen")

plt.subplot(1,2,2)
plt.imshow(np.log10(aperture_psf), cmap=colormap)
plt.xlabel("[pixels]")
plt.ylabel("[pixels]")
plt.title("Aperture PSF")

plt.figure(figsize = (14,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(input_img, cmap=colormap, extent = (-nxy*platescale/2,nxy*platescale/2,-nxy*platescale/2,nxy*platescale/2))
plt.xlabel("[arcsec]")
plt.ylabel("[arcsec]")
plt.title("Input Image")

plt.subplot(1,2,2)
plt.imshow(sensor_img, cmap=colormap)
plt.xlabel("[pixels]")
plt.ylabel("[pixels]")
plt.title("Sensor Image")

# Show center row of aperture PSF and sensor's image
plt.figure(figsize = (20,10), dpi = 100)
dx = 100 # Number of points to be plotted around the center of image
x = np.arange(center-50,center+50) # Points to be plotted

# Aperture PSF Plot
plt.subplot(1,2,1)
plt.plot(x,aperture_psf[center,x])
plt.title("Aperture PSF (Center Row)")

# Output Image Plot
plt.subplot(1,2,2)
plt.plot(x,sensor_img[center,x])
plt.title("Output Image (Center Row)")

# Show figures
plt.show()

# List distances from center of nulls in PSF
print("Aperture PSF Null Radii")
minima_i = argrelextrema(aperture_psf[center], np.less)
minima_i = np.subtract(minima_i,center)[0]
print_num = 5
minima_found = 0
for value in minima_i:
    if (value > 0) and (minima_found < print_num):
        print(value)
        minima_found = minima_found + 1
