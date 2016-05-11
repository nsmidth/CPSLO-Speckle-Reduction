# Calculates PSD with several methods and compares calculation times

# Import Modules
import sys, os
from labeyrie_classes import target, fftw_psd
import tkinter as tk
from tkinter import filedialog
from cython_psd import cython_psd

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt 
import time

# Raw data target object
raw = target()

# Start up Tkinter
root = tk.Tk()
root.withdraw()

# Prompt for FITS file locations
raw.fits.fileName = filedialog.askopenfilename(title="Select FITS file for test data")
#raw.fits.fileName = "/home/niels/Documents/FITS/KP330.fits"

# Import FITS data
raw.fits.read(numDimensions=3)
image = raw.fits.data[0]

# Number of iterations to do for benchmark
iterations = 100

# Labels for data 
time_data = {}

# Test Numpy PSD Calculation
avg_time = 0
for i in np.arange(iterations): # Repeat 
  ts = time.time()        # Start Timer
  psd = np.power(np.abs(np.fft.fft2(image)),2)   # Execute Function
  te = time.time()        # Stop timer
  avg_time += (te-ts)/iterations # Accumulate time
time_data["NumPy FFT"] = avg_time*1000

# Test SciPy PSD Calculation
avg_time = 0
for i in np.arange(iterations): # Repeat 
  ts = time.time()        # Start Timer
  psd = np.power(np.abs(sp.fftpack.fft2(image.astype(np.float32))),2)   # Execute Function
  te = time.time()        # Stop timer
  avg_time += (te-ts)/iterations # Accumulate time
time_data["SciPy FFTpack"] = avg_time*1000

# Test Optimized Scipy PSD Calculation
avg_time = 0
for i in np.arange(iterations): # Repeat 
  ts = time.time()        # Start Timer
  psd = np.abs(sp.fftpack.fft2(image.astype(np.float32)))**2   # Execute Function
  te = time.time()        # Stop timer
  avg_time += (te-ts)/iterations # Accumulate tim
time_data["Modified SciPy FFTpack"] = avg_time*1000

# Test Cython NumPy PSD Calculation
avg_time = 0
for i in np.arange(iterations): # Repeat 
  ts = time.time()        # Start Timer
  psd = cython_psd(image) # Execute function
  te = time.time()        # Stop timer
  avg_time += (te-ts)/iterations # Accumulate time
time_data["Cython NumPy FFT"] = avg_time*1000

# Test FFTW PSD Calculation
avg_time = 0
for i in np.arange(iterations): # Repeat 
  ts = time.time()        # Start Timer
  psd = fftw_psd(image)   # Execute Function
  te = time.time()        # Stop timer
  avg_time += (te-ts)/iterations # Accumulate time
time_data["Ctypes FFTW"] = avg_time*1000

# Print Data
for name, value in time_data.items():
  print('{0:25} ==> {1:10.2f}ms'.format(name,value))
  
# Graph Data
ind = np.arange(len(time_data))
width = 0.35
fig = plt.figure()
ax = fig.add_subplot(111)
ax.bar(ind, time_data.values())
ax.set_xticklabels(time_data.keys(), rotation=45)
ax.set_xticks(ind+width/2)
ax.set_ylabel("Execution Time [ms]")
ax.set_title("Execution Time Comparison of Different PSD Calculations")
ax.grid(True)
plt.tight_layout()
plt.show()
  