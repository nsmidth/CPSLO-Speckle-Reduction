import sys, os
from classes_labeyrie import target,deconvolved
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift

binary = target()
reference = target()
deconv = deconvolved()

# Define filenames for each target
binary.fits.fileName = "/home/niels/Documents/FITS/KP330.fits"
reference.fits.fileName = "/home/niels/Documents/FITS/KP331.fits"

# Import each target
binary.fits.read(numDimensions=3,printInfo=False)
reference.fits.read(numDimensions=3,printInfo=False)

# Calculate PSD of each target
binary.psdCalc()
print("Binary PSD Calc Complete")
reference.psdCalc()
print("Reference PSD Calc Complete")

## Deconvolve reference and binary stars with different methods
nxy = 512
center = int(nxy/2)
# Assign wiener variables
H = reference.psd.data
G = binary.psd.data
# Inverse filtering deconvolution
F_hat_inverse = G/H
f_hat_inverse = fftshift(np.abs(ifft2(F_hat_inverse)))

## Simplified Wiener filtering
# Set Constant
k1 = 1E-5
# Wiener filter and calculate acorr
F_hat_wiener1 = G*(1/H)*((H**2)/(H**2+k1))
f_hat_wiener1 = fftshift(np.abs(ifft2(F_hat_wiener1)))

## Wiener filtering with LPF for Signal PSD
# Setting constant and LPF radius
radius = 30
k2 = 1E-7
# Create centered meshgrid of image
xx,yy = np.meshgrid(np.arange(nxy),np.arange(nxy))
xx = np.subtract(xx,center)
yy = np.subtract(yy,center)
rr = np.power(np.power(xx,2)+np.power(yy,2),0.5)
# Create LPF filter image
lpf = np.exp(-(np.power(rr,2)/(2*np.power(radius,2))))
# Wiener filter and calculate acorr
F_hat_wiener2 = G*(1/H)*((H**2)/(H**2+k2/lpf))
f_hat_wiener2 = fftshift(np.abs(ifft2(F_hat_wiener2)))

## Display Images
colormap = "jet"

plt.figure(figsize = (16,16), dpi = 100)
plt.subplot(2,2,1)
plt.imshow(binary.fits.data[10], cmap=colormap)
plt.title("Binary Star Image")
plt.subplot(2,2,2)
plt.imshow(reference.fits.data[0], cmap=colormap)
plt.title("Reference Star Image")
plt.subplot(2,2,3)
plt.imshow(np.log10(G), cmap=colormap)
plt.title("Binary Star PSD")
plt.subplot(2,2,4)
plt.imshow(np.log10(H), cmap=colormap)
plt.title("Reference Star PSD")

plt.figure(figsize = (12,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(F_hat_wiener1, cmap=colormap)
plt.title("Deconvolved PSD")
plt.subplot(1,2,2)
plt.imshow(f_hat_wiener1, cmap=colormap)
plt.title("Deconvolved Acorr")

plt.figure(figsize = (12,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(F_hat_wiener2, cmap=colormap)
plt.title("Deconvolved PSD")
plt.subplot(1,2,2)
plt.imshow(f_hat_wiener2, cmap=colormap)
plt.title("Deconvolved Acorr")

plt.figure(figsize = (12,6), dpi = 100)
plt.subplot(1,2,1)
plt.imshow(((H**2)/(H**2+k1)), cmap=colormap)
plt.title("Simplified Wiener")
plt.subplot(1,2,2)
plt.imshow(((H**2)/(H**2+k2/lpf)), cmap=colormap)
plt.title("Wiener w/ LPF")

plt.show()