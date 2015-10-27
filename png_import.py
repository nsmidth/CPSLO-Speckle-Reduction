"""
    This is a test of importing an image

    Niels Smidth - 10/16/15
"""

# Import Modules
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

# Import images
imgSingleStar = mpimg.imread('10-15 Single Speckle Test.png')
imgDoubleStar = mpimg.imread('10-15 Binary Speckle Test.png')


# View images
imgSinglePlot = plt.imshow(imgSingleStar)
imgDoublePlot = plt.imshow(imgDoubleStar)

plt.show()

# Debug End of Code
print("Done Running")
