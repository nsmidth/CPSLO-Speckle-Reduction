#
#
# Testing some principles of normalization
#
#

import numpy as np
from scipy.fftpack import fft2, ifft2

# Blank array, take FFT of it
a = np.array([ [1,1,1], [1,1,1], [1,1,1] ])
A = fft2(a)

print("a = ","\n",a)
print("A = ","\n",A)
print("Needs a factor of normalization of 1/(elements in a)")
print()

# Showing normalization
A_norm = fft2(a)/a.size

print("Normalized A = A / (A.size) = ", "\n", A_norm)
print()

# PSD of normalized and non-normalized
a_PSD = np.abs(A)**2
a_PSD_norm = np.abs(A_norm)**2

print("PSD of (a not normalized) = ", "\n", a_PSD)
print("PSD of (a normalized) = ", "\n", a_PSD_norm)
print()

# Autocorrelation of normalized and non-normalized
a_autocorr = ifft2(a_PSD)
a_autocorr_norm = ifft2(a_PSD_norm)

a_autocorr_norm2 = ifft2(a_PSD_norm)*a_PSD_norm.size
a_autocorr_norm3 = ifft2(a_PSD_norm)*(a_PSD_norm.size**2)

print("Autocorrelation of (a not normalized) = ", "\n", a_autocorr)
print("Autocorrelation of (a normalized) = ", "\n", a_autocorr_norm)
print("We see that the ifft function performs the normalization function again")
print(" Because we applied normalization before the squaring of PSD, this function")
print(" was squared during PSD calculation. So we need to multiply by ")
print(" a_PSD_norm.size^2 in order to get autocorr back on scale ")
print(" with original value, or just .size if we want the autocorr ")
print(" to be normalized to 1")

print("Autocorrelation of (a normalized) * size = ", "\n", a_autocorr_norm2)
print("Autocorrelation of (a normalized) * size^2 = ", "\n", a_autocorr_norm3)
print()
