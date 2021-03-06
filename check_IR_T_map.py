
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from basic_functions import *
from cmb import *
from flat_map import *

import sys


pathFig = "./figures/sehgal_maps/ir_sources/"

# read and add up all my IR temperature maps [muK]
nSide = 4096
tMap = np.zeros(hp.nside2npix(nSide))
#SourceCatalogs = ['IRgal_S_'+str(i) for i in range(1, 12)] + ['IRBlastPop']
SourceCatalogs = ['IRgal_S_'+str(i) for i in range(1, 11)] + ['IRBlastPop']

for sourceCatalog in SourceCatalogs:
   # read the healpy map
   path = "./output/sehgal_maps/ir_sources/t_ir_"+sourceCatalog+"_sehgal_148ghz_muk.fits"
   tMap += hp.read_map(path)

# read the official Sehgal map [mJy/sr]
path = "./input/sehgal_maps/148_ir_pts_healpix.fits"
sehgalTMap = hp.read_map(path)
sehgalTMap = hp.ud_grade(sehgalTMap, 4096)
# convert from [mJy /sr] to [dT/Tcmb]
sehgalTMap /= 1.072480e9
# convert to muKcmb
sehgalTMap *= 2.726e6


# Plot my t map
hp.mollview(tMap, title='My IR T map')
plt.savefig(pathFig+"mytmap.pdf")
plt.show()

# Plot the official Sehgal t map
hp.mollview(sehgalTMap, title='Sehgal IR T map')
plt.savefig(pathFig+"sehgaltmap.pdf")
plt.show()

# compare the maps
diff = np.log10(np.abs(tMap-sehgalTMap))
hp.mollview(diff, title='My IR T map VS Sehgal T map')
plt.savefig(pathFig+"logabsdiff_mytmap_vs_sehgaltmap.pdf")
plt.show()

