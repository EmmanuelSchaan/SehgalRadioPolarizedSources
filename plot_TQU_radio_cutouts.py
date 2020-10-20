import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from basic_functions import *
from cmb import *
from flat_map import *

# needed on lawrencium
plt.switch_backend('Agg')

####################################################################
# Choose frequency

nu = np.float(sys.argv[1])
nuStr = str(np.int(nu/1.e9))
fluxCutmJy = np.float(sys.argv[2])
lKnee = np.float(sys.argv[3])
aKnee = np.float(sys.argv[4])
beamFwhm = np.float(sys.argv[5])
noiseT = np.float(sys.argv[6])

print "nu", nu
print "nuStr", nuStr
print "fluxCutmJy", fluxCutmJy
print "lKnee", lKnee
print "aKnee", aKnee
print "beamFwhm", beamFwhm
print "noiseT", noiseT


##nu = 90.e9  # [Hz]
##nuStr = '90'
#nu = 148.e9  # [Hz]
#nuStr = '148'
#
#
#fluxCutmJy = 10.  #5.   #2.  # [mJy]
#
#
## wide survey
#lKnee = 700.
#aKnee = 1.4
#beamFwhm = 1.4 # [arcmin] at 148 GHz
##beamFwhm = 2.2 # [arcmin] at 90 GHz
#noiseT = 2.  # [muK*arcmin] at 148 GHz 
##noiseT = 2.  # [muK*arcmin] at 90 GHz
#
### deep survey
##lKnee = 200.
##aKnee = 2.
##beamFwhm = 1.5 # [arcmin] at 148 GHz
###beamFwhm = 2.3 # [arcmin] at 90 GHz
##noiseT = 0.96  # [muK*arcmin] at 148 GHz 
###noiseT = 0.68  # [muK*arcmin] at 90 GHz


cmb = CMB(beam=beamFwhm, noise=noiseT, nu1=nu, nu2=nu, lMin=30., lMaxT=3.e3, lMaxP=5.e3, fg=True, atm=True, atmProp=[lKnee, aKnee, lKnee, aKnee], name=None)


####################################################################
# # Extract as many square cutouts as possible

# cutout dimensions
# map side in lon and lat
dLon = 10.# [deg]
dLat = 10.# [deg]
lonRange = np.array([-dLon/2., dLon/2.]) # [deg]
latRange = np.array([-dLat/2., dLat/2.]) # [deg]
pixRes = 0.5/60.  #0.5 / 60.  # [arcmin] to [deg]
# number of pixels on the side
xSize = np.int(np.ceil(dLon / pixRes))
ySize = np.int(np.ceil(dLat / pixRes))

baseMap = FlatMap(nX=xSize, nY=ySize, sizeX=dLon*np.pi/180., sizeY=dLat*np.pi/180.)


####################################################################

nPatches = 0

pathOut = "./output/sehgal_maps/radio_sources/cutouts/"
pathOut + "ps_official_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt", cutSehgalTMap
pathOut + "ps_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt", cutTMap
pathOut + "ps_sehgal_"+nuStr+"ghz_Q_patch"+str(nPatches)+".txt", cutQMap
pathOut + "ps_sehgal_"+nuStr+"ghz_U_patch"+str(nPatches)+".txt", cutUMap
pathOut + "kappa_sehgal_patch"+str(nPatches)+".txt", cutKappaMap
pathOut + "ps_mask_"+nuStr+"ghz_"+str(np.int(round(fluxCutmJy)))+"mJy_beam"+str(round(beamFwhm,1))+"_noise"+str(round(noiseT,2))+"_lknee"+str(np.int(lKnee))+"_aknee"+str(round(aKnee,1))+"_T_patch"+str(nPatches)+".txt", psMask


