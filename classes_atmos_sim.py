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

class atmospheric_simulation():
  def __init__(self):
    # Aperture/Telescope Specifications:
    self.diameter_m = 2.133 # Mirror diameter in meters
    self.focal_length = 129.69 # Effective focal length in meters

    # Camera Specs
    self.wavelength = 0.8E-6 # Wavelength of light
    self.pixel = 8E-6 # length of pixel side in m
    self.bits = 14 # Bits in camera ADC
    self.nxy = 512 # Length of sensor side in pixels
    self.center = int(nxy/2) # Center of sensor in pixels
    self.platescale = 206265*self.pixel/(self.focal_length) # Calculate plate scale
    self.gamma = 1.6 # Gamma correction

    # Binary star specs
    self.rho = 1.5 # Set separation in arcseconds
    self.phi = 45 # Set angle in degrees  
    
    # Atmospheric Specs
    self.alpha = 1/100 # Multiplicative constant
    self.r0 = 0.2 # Fried Parameter
    
  def create_input_image(self): 
    # Calculate coordinates of stars
    x = int( self.rho/(2*self.platescale) * np.cos(np.deg2rad(self.phi)) )
    y = int( self.rho/(2*self.platescale) * np.sin(np.deg2rad(self.phi)) )
    x1 = self.center + x
    y1 = self.center + y
    x2 = self.center - x
    y2 = self.center - y    
    # Empty input image
    self.input_img = np.zeros((nxy,nxy)) 
    # Place stars on image
    self.input_img[y1,x1] = 1 
    self.input_img[y2,x2] = 1
    # Scale image power to 1
    input_img_power = np.sum(np.power(self.input_img,2))
    self.input_img = np.divide(self.input_img,np.sqrt(input_img_power))
    
  def get_psf(self):

    ## Telescope aperture creation:
    # Total spatial sample range
    X_aperture_s = 1/self.pixel 
    # dx value for sampled aperture image
    dx_aperture_s = X_aperture_s/nxy 
    # Coordinates of sampled image
    x_aperture_s = np.arange(0,X_aperture_s,dx_aperture_s) - X_aperture_s/2
    # Meshgrid of sampled coordinates
    xx_s,yy_s = np.meshgrid(x_aperture_s,x_aperture_s)
    # Scaled aperture diameter to effectively resample aperture image
    diameter_s = self.diameter_m/(self.focal_length*self.wavelength)
    # Draw new circle at correct dimensions
    # Calculate grid of distances from circle center
    circle_s = (xx_s) ** 2 + (yy_s) ** 2 
    # Draw boolean circle
    circle_s = circle_s < (diameter_s/2)**2 
    # Convert boolean circle to int
    circle_s= circle_s.astype(np.int64)
    # Scale aperture image power to 1
    aperture_screen_power = np.sum(np.power(circle_s,2))
    # Save aperture image in units of meters
    self.aperture_screen_s = np.divide(circle_s,np.sqrt(aperture_screen_power))
    # Calculate effective size of sampled aperture image in meters
    X_aperture_s_meff = self.focal_length*self.wavelength/self.pixel

    ## Phase screen creation:
    # Generate random image
    # To be used in creating random atmospheric element
    phase_phase = np.multiply(np.pi,np.random.normal(loc=0,scale=1,size=(nxy,nxy)))
    # Total array sample size
    d_aperture = X_aperture_s_meff
    # Spatial sample resolution
    dxy = d_aperture/nxy
    # Spatial frequency resolution
    df = 1/(d_aperture) 
    # Image sample indices array
    x = np.multiply( np.subtract(np.arange(self.nxy),self.center), dxy ) 
    # Spatial Frequency indices array
    xf = np.multiply( np.subtract(np.arange(self.nxy),self.center), df ) 
    # Meshgrid of spatial frequency domain
    [xx,yy]=np.meshgrid(xf,xf)
    # Radius from center meshgrid
    rr = (np.sqrt(np.power(xx,2)+np.power(yy,2)))
    # Calculate Kolmogorov spectral density
    phase_PSD = np.power(rr,-11/3)
    phase_PSD = np.multiply(self.alpha*0.023/(self.r0**(5/3)),phase_PSD)
    # Set DC component to 0 (previous calc attempts to set to 1/0)
    phase_PSD[self.center,self.center] = 0 
    # Construct phase screen spectrum
    phase_screen_f = np.multiply(np.sqrt(phase_PSD),np.exp(1j*phase_phase))
    # Calculate phase screen
    phase_screen = np.real(ifft2(fftshift(phase_screen_f)*nxy*nxy))
    # Create complex atmospheric screen
    self.atmosphere_screen = np.exp(np.multiply(1j,phase_screen))

    # Generate total screen, combining atmosphere and aperture
    self.pupil_screen = np.multiply(self.atmosphere_screen,self.aperture_screen_s)

    ## Calculate system's total response 
    # Calculate total PSF of system
    self.psf = fftshift(fft2(self.pupil_screen))
    self.psf = np.power(np.abs(self.psf),2)
    # Normalize PSF
    psf_power = np.sum(np.power(self.psf,2))
    self.psf = np.divide(self.psf,psf_power)
    # Gamma correct PSF
    self.psf = np.power(self.psf,1/self.gamma)  
                     
  def get_binary(self):
        # Convolve PSF with input image using FFT
    sensor_img = fftconvolve(input_img,psf)
    # Save the center 512x512 image
    sensor_img = sensor_img[center:center+nxy,center:center+nxy]
  def add_noise(self,img):
    pass
    
    

  
  
  ## Adding noise to images
def add_gaus_noise(img, var):
    # Size of input image?
    nxy = np.shape(img)[0]
    
    # Array to be returned
    img_noisy = np.zeros((nxy,nxy))    
    
    # Create noise image
    noise = np.random.normal(loc=0,scale=var,size=(nxy,nxy))
    # Add noise image to simulated image
    img_noisy = img+noise
    # Turn all negative pixels to 0
    img_noisy[img_noisy<0] = 0
    
    # Return noisy image
    return img_noisy
  
def add_shot_noise(img, photons):
    # Size of input image?
    nxy = np.shape(img)[0]    
    
    # Empty arrays to use
    img_normalized = np.zeros((nxy,nxy))
    img_scaled = np.zeros((nxy,nxy))
    img_shot_noise = np.zeros((nxy,nxy))
    
    # Normalize input image
    img_normalized = np.abs(img)/(np.sum(img))
    # Scale image to number of photons desired
    img_scaled = img_normalized*photons
    # Calculate image with shot noise
    img_shot_noise = np.random.poisson(lam=img_scaled, size=None)

    # Return noisy image
    return img_shot_noise