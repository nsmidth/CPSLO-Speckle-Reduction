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
from speckle_fns import deconv0

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
zeroSub = 0.001
psdDeconv0 = deconv0(psdDoubleStar, psdSingleStar, zeroSub)

# Perform postprocessing on images
acorrSingleStar = postprocess( psdSingleStar )
acorrDoubleStar = postprocess( psdDoubleStar )
#acorrDeconv0 = postprocess( psdDeconv0 )

# View images
printImg = input("Print star images (y/n)?   ")
if (printImg.lower() == "y"):
    print_star(imgSingleStar, psdSingleStar, acorrSingleStar)
    print_star(imgDoubleStar, psdDoubleStar, acorrDoubleStar)
    #print_star(imgDoubleStar, psdDeconv0, acorrDeconv0)

# Debug End of Code
print("Done Running")
