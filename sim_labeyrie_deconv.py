# Includes
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.signal import fftconvolve
from scipy.signal import argrelextrema
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

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
n_exposures = 100

## Creating Binary Star Input Image
# For KP 2.1m telescope, input image is in units 
#  of 25.44 milliarcseconds
platescale = 0.02544 # Plate scale in arcsec/pixel
rho = 1.5 # Set separation in arcseconds
phi = 45 # Set angle in degrees
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
x_aperture_s = np.arange(0,X_aperture_s,dx_aperture_s) - X_aperture_s/2
# Meshgrid of sampled coordinates
xx_s,yy_s = np.meshgrid(x_aperture_s,x_aperture_s)
# Scaled aperture diameter to 
#  effectively resample aperture image
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

# Create image accumulation array
psf_avg = np.zeros((nxy,nxy))
sensor_img_psd_avg = np.zeros((nxy,nxy))
atmosphere_psd_avg = np.zeros((nxy,nxy))

# Generate multiple images, integrate them
for i in np.arange(n_exposures):
    # Generate random phase
    phase_phase = np.multiply(np.pi,np.random.normal(loc=0,scale=1,size=(nxy,nxy)))
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
    # Integrate PSFs
    psf_avg += psf
    
    # Calculate PSD of binary star image
    sensor_img_psd = np.power(np.abs(fftshift(fft2(sensor_img))),2)
    # Integrate PSDs
    sensor_img_psd_avg += sensor_img_psd
    
    # Calculate PSD of PSF
    atmosphere_psd = np.power(np.abs(fftshift(fft2(psf))),2)
    # Integrate PSDs
    atmosphere_psd_avg += atmosphere_psd
    
# Calculate average of PSF/PSDs
psf_avg /= n_exposures
sensor_img_psd_avg /= n_exposures
atmosphere_psd_avg /= n_exposures

# Calculate PSDs
input_img_psd = np.power(np.abs(fftshift(fft2(input_img))),2)

# Calculate Acorrs
input_img_acorr = np.abs(fftshift(ifft2(input_img_psd)))
atmosphere_acorr_avg = np.abs(fftshift(ifft2(atmosphere_psd_avg)))
sensor_img_acorr_avg = np.abs(fftshift(ifft2(sensor_img_psd_avg)))

# Convolve small circles with input images to make them more visible
# Create meshgrid
xx,yy = np.meshgrid(np.arange(nxy),np.arange(nxy))
# Calculate grid of distances from circle center
radius = (xx-center) ** 2 + (yy-center) ** 2 
# Draw boolean circle
circle = (radius < (3)**2).astype(np.int64)
# Convolve circle with difficult to see images
input_img_emph = fftconvolve(input_img,circle)[center:center+nxy,center:center+nxy]
input_img_acorr_emph = fftconvolve(input_img_acorr,circle)[center:center+nxy,center:center+nxy]

# Deconvolve reference from binary
deconvolved_psd_avg = np.divide(sensor_img_psd_avg, atmosphere_psd_avg)
# LPF Filter
lpf_radius = 70
lpf = np.exp(-(np.power(radius,2)/(2*np.power(lpf_radius,4))))
#deconvolved_psd_filtered_avg = lpf*deconvolved_psd_avg

deconvolved_acorr_avg = np.abs(fftshift(ifft2(deconvolved_psd_avg)))
#deconvolved_acorr_filered_avg = np.abs(fftshift(ifft2(deconvolved_psd_filtered_avg)))

colormap = "jet"

plt.figure(figsize = (18,8), dpi = 150)
plt.subplot(1,2,1)
plt.imshow(psf, cmap=colormap)
plt.title("Short Exposure PSF")
plt.subplot(1,2,2)
plt.imshow(psf_avg, cmap=colormap)
plt.title("Long Exposure PSF")

plt.figure(figsize = (16,16), dpi = 50)
plt.subplot(3,3,1)
plt.imshow(input_img_emph, cmap=colormap)
plt.title("Binary Star Object")
plt.subplot(3,3,2)
plt.imshow(psf, cmap=colormap)
plt.title("Atmospheric/Aperture PSF")
plt.subplot(3,3,3)
plt.imshow(sensor_img, cmap=colormap)
plt.title("Binary Star Image")
plt.subplot(3,3,4)
plt.imshow(np.log10(input_img_psd), cmap=colormap)
plt.title("Binary Star Object PSD")
plt.subplot(3,3,5)
plt.imshow(np.log10(atmosphere_psd_avg), cmap=colormap)
plt.title("Avg PSF PSD")
plt.subplot(3,3,6)
plt.imshow(np.log10(sensor_img_psd_avg), cmap=colormap)
plt.title("Avg Binary Star Image PSD")
plt.subplot(3,3,7)
plt.imshow(input_img_acorr_emph, cmap=colormap)
plt.title("Binary Star Object Autocorrelation")
plt.subplot(3,3,8)
plt.imshow(atmosphere_acorr_avg, cmap=colormap)
plt.title("Avg PSF Autocorrelation")
plt.subplot(3,3,9)
plt.imshow(sensor_img_acorr_avg, cmap=colormap)
plt.title("Avg Binary Star Image Autocorrelation")

fig = plt.figure(figsize = (10,10), dpi = 100)
ax = fig.gca(projection='3d')
[xx,yy] = np.meshgrid(np.arange(nxy), np.arange(nxy))
ax.plot_surface(X=xx,
                Y=yy,
                Z=np.log10(atmosphere_psd_avg),
                cmap="jet",
                linewidth=0, 
                antialiased=True)
plt.title("Reference Star PSD")

fig = plt.figure(figsize = (16,6), dpi = 100)
ax = fig.add_subplot(1,2, 1, projection='3d')
ax.plot_wireframe(  X=xx,
                    Y=yy,
                    Z=np.log10(1/atmosphere_psd_avg+0.001),
                    cmap="jet",
                    linewidth=0.1, 
                    antialiased=True)
ax.set_title('1/(Reference Star PSD)')
ax = fig.add_subplot(1,2, 2, projection='3d')
ax.plot_wireframe(  X=xx,
                    Y=yy,
                    Z=np.log10(lpf),
                    cmap="jet",
                    linewidth=0.1, 
                    antialiased=True)
ax.set_title('Lowpass Filter')

fig = plt.figure(figsize = (8,8), dpi = 100)
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_wireframe(  X=xx,
                    Y=yy,
                    Z=np.log10(lpf/atmosphere_psd_avg+0.001),
                    cmap="jet",
                    linewidth=0.1, 
                    antialiased=True)
ax.set_title('Lowpass Filter/(Reference Star PSD)')

plt.figure(figsize = (14,14), dpi = 100)
plt.subplot(2,2,1)
plt.imshow(np.log10(sensor_img_psd_avg), cmap=colormap)
plt.title("Avg Binary Star Image PSD")
plt.subplot(2,2,2)
plt.imshow( np.log10(deconvolved_psd_avg), cmap=colormap)
plt.title("Avg Deconvolved PSD")
plt.subplot(2,2,3)
plt.imshow(sensor_img_acorr_avg, cmap=colormap)
plt.title("Avg Binary Star Image Autocorrelation")
plt.subplot(2,2,4)
plt.imshow(deconvolved_acorr_avg, cmap=colormap)
plt.title("Avg Deconvolved Autocorrelation")

import sys, os
from classes_labeyrie import target,deconvolved
import tkinter as tk
from tkinter import filedialog

binary = target()
reference = target()
deconv = deconvolved()

# Define filenames for each target
binary.fits.fileName = "/home/niels/Documents/FITS/KP336.fits"
reference.fits.fileName = "/home/niels/Documents/FITS/KP338.fits"

# Import each target
binary.fits.read(numDimensions=3,printInfo=False)
reference.fits.read(numDimensions=3,printInfo=False)

# Calculate PSD of each target
binary.psdCalc()
print("Binary PSD Calc Complete")
reference.psdCalc()
print("Reference PSD Calc Complete")

# Deconvolve reference and binary stars
deconv.psdDeconvolve(binary.psd.data,reference.psd.data,1e-12)

# Perform filtering on output star object
deconv.psdFilter(lpfRadius = 20, interference=False)

# Get autocorrelogram
deconv.acorrCalc()

plt.figure(figsize = (12,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(binary.fits.data[10], cmap=colormap)
plt.title("Binary Star Image")
plt.subplot(1,2,2)
plt.imshow(reference.fits.data[0], cmap=colormap)
plt.title("Reference Star Image")

plt.figure(figsize = (6,12), dpi = 100)
plt.subplot(2,1,1)
plt.imshow(deconv.psdFiltered.data, cmap=colormap)
plt.title("Filtered Deconvolved PSD")
plt.subplot(2,1,2)
plt.imshow(deconv.acorr.data, cmap=colormap)
plt.title("Autocorrelation")

plt.show()
