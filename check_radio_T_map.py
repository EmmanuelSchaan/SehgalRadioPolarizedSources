
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from basic_functions import *
from cmb import *
from flat_map import *

import sys


pathFig = "./figures/sehgal_maps/radio_sources/"


# read the healpy map
path = "./output/sehgal_maps/radio_sources/t_radio_sehgal_148ghz_muk.fits"
tMap = hp.read_map(path)

# read the official Sehgal map [mJy/sr]
path = "./input/sehgal_maps/148_rad_pts_healpix.fits"
sehgalTMap = hp.read_map(path)
sehgalTMap = hp.ud_grade(sehgalTMap, 4096)
# convert from [mJy /sr] to [dT/Tcmb]
sehgalTMap /= 1.072480e9
# convert to muKcmb
sehgalTMap *= 2.726e6


# Plot my t map
hp.mollview(np.log10(tMap), title='My radio T map')
plt.savefig(pathFig+"mytmap.pdf")
plt.show()

# Plot the official Sehgal t map
hp.mollview(np.log10(sehgalTMap), title='Sehgal radio T map')
plt.savefig(pathFig+"sehgaltmap.pdf")
plt.show()

# compare the maps
diff = np.log10(np.abs(tMap-sehgalTMap))
hp.mollview(diff, title='My radio T map VS Sehgal T map')
plt.savefig(pathFig+"logabsdiff_mytmap_vs_sehgaltmap.pdf")
plt.show()

