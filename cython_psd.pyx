# Cython module for calculating PSD
cimport numpy as np
import numpy as np

def cython_psd(np.ndarray image):

  cdef np.ndarray psd = np.zeros([512,512], dtype=np.float32)

  # Calculate 2D power spectrum
  # This gives us only real values
  psd = np.abs(np.fft.fft2(image.astype(np.float32)))**2
  
  return psd
