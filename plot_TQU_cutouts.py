import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from basic_functions import *
from cmb import *

import flat_map
reload(flat_map)
from flat_map import *

# needed on lawrencium
# not on cori
#plt.switch_backend('Agg')

####################################################################
# Choose frequency

#nu = np.float(sys.argv[1])
#nuStr = str(np.int(nu/1.e9))
#fluxCutmJy = np.float(sys.argv[2])
#lKnee = np.float(sys.argv[3])
#aKnee = np.float(sys.argv[4])
#beamFwhm = np.float(sys.argv[5])
#noiseT = np.float(sys.argv[6])
#
#print "nu", nu
#print "nuStr", nuStr
#print "fluxCutmJy", fluxCutmJy
#print "lKnee", lKnee
#print "aKnee", aKnee
#print "beamFwhm", beamFwhm
#print "noiseT", noiseT

# Foreground map frequency
#
nu = 90.e9  # [Hz]
nuStr = '90'
#nu = 148.e9  # [Hz]
#nuStr = '148'



# Mask properties
#
fluxCutmJy = 2. #10.  #5.   #2.  # [mJy]

# wide survey
lKnee = 700.
aKnee = 1.4
#beamFwhm = 1.4 # [arcmin] at 148 GHz
beamFwhm = 2.2 # [arcmin] at 90 GHz
#noiseT = 2.  # [muK*arcmin] at 148 GHz 
noiseT = 2.  # [muK*arcmin] at 90 GHz

## deep survey
#lKnee = 200.
#aKnee = 2.
#beamFwhm = 1.5 # [arcmin] at 148 GHz
##beamFwhm = 2.3 # [arcmin] at 90 GHz
#noiseT = 0.96  # [muK*arcmin] at 148 GHz 
##noiseT = 0.68  # [muK*arcmin] at 90 GHz


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
####################################################################
# read one patch

nPatches = 1


####################################################################
# radio maps

pathOutRadio = "./output/sehgal_maps/radio_sources/cutouts/"
cutRadioSehgalTMap = np.genfromtxt(pathOutRadio + "ps_official_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt")
cutRadioTMap = np.genfromtxt(pathOutRadio + "ps_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt")
cutRadioQMap = np.genfromtxt(pathOutRadio + "ps_sehgal_"+nuStr+"ghz_Q_patch"+str(nPatches)+".txt")
cutRadioUMap = np.genfromtxt(pathOutRadio + "ps_sehgal_"+nuStr+"ghz_U_patch"+str(nPatches)+".txt")


psRadioMask = np.genfromtxt(pathOutRadio + "ps_mask_"+nuStr+"ghz_"+str(np.int(round(fluxCutmJy)))+"mJy_beam"+str(round(beamFwhm,1))+"_noise"+str(round(noiseT,2))+"_lknee"+str(np.int(lKnee))+"_aknee"+str(round(aKnee,1))+"_T_patch"+str(nPatches)+".txt")

####################################################################

cutKappaMap = np.genfromtxt(pathOutRadio + "kappa_sehgal_patch"+str(nPatches)+".txt")

####################################################################
# IR maps
'''
pathOutIR = "./output/sehgal_maps/ir_sources/cutouts/"

# read cutouts for all IR source catalogs
cutIRTMap = np.zeros((10, baseMap.nX, baseMap.nY))
cutIRQMap = np.zeros((10, baseMap.nX, baseMap.nY))
cutIRUMap = np.zeros((10, baseMap.nX, baseMap.nY))
for i in range(10):
   sourceCatalog = "IRgal_S_"+str(i+1)

   #cutIRSehgalTMap[i,:,:] = np.genfromtxt(pathOutIR + "ir_official_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt")
   cutIRTMap[i,:,:] = np.genfromtxt(pathOutIR + "ir_"+sourceCatalog+"_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt")
   cutIRQMap[i,:,:] = np.genfromtxt(pathOutIR + "ir_"+sourceCatalog+"_sehgal_"+nuStr+"ghz_Q_patch"+str(nPatches)+".txt")
   cutIRUMap[i,:,:] = np.genfromtxt(pathOutIR + "ir_"+sourceCatalog+"_sehgal_"+nuStr+"ghz_U_patch"+str(nPatches)+".txt")

   #psIRMask[i,:,:] = np.genfromtxt(pathOutIR + "ir_"+sourceCatalog+"_mask_"+nuStr+"ghz_"+str(np.int(round(fluxCutmJy)))+"mJy_beam"+str(round(beamFwhm,1))+"_noise"+str(round(noiseT,2))+"_lknee"+str(np.int(lKnee))+"_aknee"+str(round(aKnee,1))+"_T_patch"+str(nPatches)+".txt")
'''

####################################################################


pathFig = "./figures/summary_figures/"


def plotMap(d, title=None, unit=None, save=False, path=None):
   dF = baseMap.fourier(d)
   f = lambda l: baseMap.gaussianBeam(l, 2.*np.pi/180./60.)
   dF = baseMap.filterFourierIsotropic(f, dF)
   d = baseMap.inverseFourier(dF)
   s = np.std(d)
   baseMap.plot(d, vlim=(-s, s), cmap='bwr', unit=unit, title=title, save=save, path=path)


path = pathFig + "ps_sehgal_"+nuStr+"ghz_Q_patch"+str(nPatches)+".pdf"
plotMap(cutRadioQMap.copy(), title='Radio '+nuStr+' GHz', unit=r'$\mu$K', save=True, path=path)
'''
path = pathFig + "ir_"+"IRgal_S_1-10"+"_sehgal_"+nuStr+"ghz_Q_patch"+str(nPatches)+".pdf"
plotMap(np.sum(cutIRQMap, axis=0), title='IR '+nuStr+' GHz', unit=r'$\mu$K', save=True, path=path)
'''
path = pathFig + "kappa_sehgal_patch"+str(nPatches)+".pdf"
plotMap(cutKappaMap.copy(), title=r'$\kappa_\text{CMB}$', save=True, path=path)
