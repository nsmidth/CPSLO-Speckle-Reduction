"""
    To test deconvolution, it is easiest to use some image files
    that I know will produce nice deconvolution results. 

    Niels Smidth - 10/16/15
"""

# Import Modules
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from speckle_fns import preprocess,postprocess,print_star
from speckle_fns import deconv1

# Import images
imgSingleStar = mpimg.imread('10-15 Single Speckle Test.png')
imgDoubleStar = mpimg.imread('10-15 Binary Speckle Test.png')

# Only use one channel of these B/W images
imgSingleStar = imgSingleStar[:,:,1]
imgDoubleStar = imgDoubleStar[:,:,1]

# Perform preprocessing on both images
psdSingleStar = preprocess( imgSingleStar )
psdDoubleStar = preprocess( imgDoubleStar )

# Perform deconvolution on images
constant = 0.001
psdDeconv1 = deconv1(psdDoubleStar, psdSingleStar, constant)

# Perform postprocessing on images
acorrSingleStar = postprocess( psdSingleStar )
acorrDoubleStar = postprocess( psdDoubleStar )
acorrDeconv1 = postprocess( psdDeconv1 )

# View images
printImg = input("Print star images (y/n)?   ")
if (printImg.lower() == "y"):
    print_star(imgSingleStar, psdSingleStar, acorrSingleStar)
    print_star(imgDoubleStar, psdDoubleStar, acorrDoubleStar)
    print_star(imgDoubleStar, psdDeconv1, acorrDeconv1)

# Debug End of Code
print("Done Running")
