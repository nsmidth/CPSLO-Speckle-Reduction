"""
    Experimenting with freq vs spatial domain transforms
"""

# Import Modules
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from speckle_fns import preprocess,postprocess,print_star
from speckle_fns import deconv0
from scipy.fftpack import fftshift


# Generate empty frequency spectra
freq = np.zeros( ( 32,32 ) )
freq[:,:] = 0.001

space = postprocess(freq)

plt.subplot(1,3,1)
plt.title("Freq")
plt.imshow(freq)
plt.grid()

plt.subplot(1,3,2)
plt.title("Freq Shifted")
plt.imshow(fftshift(freq))
plt.grid()

plt.subplot(1,3,3)
plt.title("Space")
plt.imshow(space)
plt.colorbar()
plt.show()
