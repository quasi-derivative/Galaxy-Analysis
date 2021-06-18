# MODULE IMPORTS
import numpy as np
from scipy.stats import norm
from astropy.wcs import WCS
from astropy.io import fits, ascii
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.convolution import Gaussian2DKernel, convolve
from reproject import reproject_interp

from mpdaf.obj import Cube, Image
from mpdaf.obj import iter_ima, iter_spe

from photutils import CircularAperture, SkyCircularAperture, aperture_photometry
import matplotlib.pyplot as plt

# set variable c to the coordinates of the galaxy (currently NGC0628)
c = SkyCoord('01h36m41.7s','+15d47m01s', frame='icrs')
# NGC0628: '01h36m41.7s','+15d47m01s'
# NGC4579: '12h37m43.5220s', '+11d49m05.498s'

# import the near infrared (7625 Å) filter image
image_i = Image('http://dustpedia.astro.noa.gr/Data/GetImage?imageName=NGC0628_SDSS_i.fits&instrument=SDSS')
image_i = image_i.subimage(center=(c.dec.deg,c.ra.deg), size=600.0)
print(image_i.shape)

# import the green (4770 Å) filter image
image_g = Image('http://dustpedia.astro.noa.gr/Data/GetImage?imageName=NGC0628_SDSS_g.fits&instrument=SDSS')
image_g = image_g.subimage(center=(c.dec.deg,c.ra.deg), size=600.0)
print(image_g.shape)

# find ratio between the two, and save into an mpdaf.obj.Image (very much hacky)
diff = image_g.data / image_i.data
imavg = image_i
imavg.data = diff

# f = plt.figure(figsize=(6,6))
# f.add_subplot(1,2, 1)
# plots the image
imavg.plot(scale='arcsinh', cmap='inferno', colorbar='v',vmin=0.24,vmax=0.56)

# This stuff is to try and do gaussian fits, you can ignore
# calculate data from the image
print(imavg.background())
print(image_g.background())
peak = imavg.peak(center=(c.dec.deg,c.ra.deg))
fwhm = imavg.fwhm(center=(c.dec.deg,c.ra.deg))
print(fwhm)
gfit = imavg.gauss_fit(maxiter=150, plot=False, center=(c.dec.deg,c.ra.deg))
print(gfit.fwhm)
print(gfit.peak)

# f.add_subplot(1,2, 2)
a = gfit.peak
b = ((gfit.center[0]+gfit.center[1]) / 2)
sigma = ((gfit.fwhm[0]+gfit.fwhm[1]) / 2) / 2.35482004503

x = np.linspace(0,1333,1333)
y = a * np.exp((x-b)**2/2*sigma**2)
# plt.plot(x, y)

plt.show()
