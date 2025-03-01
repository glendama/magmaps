import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian2DKernel
from scipy import ndimage

# Reading and showing the original map
data = np.fromfile('map1.bin', '<f4')
# Determining image size NxN
N = int(np.sqrt(len(data))) 
print (N)

image = data.reshape((N, N))
# Tuning appropriate values of vmin and vmax
plt.imshow (image, origin="lower", vmin = image.mean()*0.25, vmax = image.mean()*4)
plt.show()

## Creating a 2D Gaussian with a radius of 100 pixels  
gauss = Gaussian2DKernel(100)

## Performing the convolution and saving the results
cimage = convolve_fft(image, gauss, allow_huge=True)
np.save('map1conv', cimage)
plt.imshow(cimage, origin="lower", vmin = image.mean()*0.3, vmax = image.mean()*3)
plt.colorbar()
plt.savefig('map1conv.png')
plt.show()
