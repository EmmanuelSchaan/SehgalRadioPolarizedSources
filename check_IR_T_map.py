
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from basic_functions import *
from cmb import *
from flat_map import *

import sys


pathFig = "./figures/"

# read and add up all my IR temperature maps [muK]
nSide = 4096
tMap = np.zeros(hp.nside2npix(nSide))
SourceCatalogs = ['IRgal_S_'+str(i) for i in range(1, 12)] + ['IRBlastPop']

for sourceCatalog in SourceCatalogs:
   # read the healpy map
   path = "./output/sehgal_maps/ir_sources/t_ir_"+sourceCatalog+"_sehgal_148ghz_muk.fits"
   tMap += hp.read_map(path)

# read the official Sehgal map [mJy]
path = "./input/sehgal_maps/148_ir_pts_healpix.fits"
sehgalTMap = hp.read_map(path)


hp.mollweide(tMap-sehgalTMap, title='My IR T map VS Sehgal T map')
plt.savefig()
plt.show()

