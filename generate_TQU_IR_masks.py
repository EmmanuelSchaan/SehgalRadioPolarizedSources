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

sourceCatalog = sys.argv[1]
nu = np.float(sys.argv[2])
nuStr = str(np.int(nu/1.e9))
fluxCutmJy = np.float(sys.argv[3])
lKnee = np.float(sys.argv[4])
aKnee = np.float(sys.argv[5])
beamFwhm = np.float(sys.argv[6])
noiseT = np.float(sys.argv[7])


print("Analyzing catalog "+sourceCatalog)

print("nu", nu)
print("nuStr", nuStr)
print("fluxCutmJy", fluxCutmJy)
print("lKnee", lKnee)
print("aKnee", aKnee)
print("beamFwhm", beamFwhm)
print("noiseT", noiseT)


#####################################################################
# CMB power spectra

cmb = CMB(beam=beamFwhm, noise=noiseT, nu1=nu, nu2=nu, lMin=30., lMaxT=3.e3, lMaxP=5.e3, fg=True, atm=True, atmProp=[lKnee, aKnee, lKnee, aKnee], name=None)

####################################################################
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
# Generate and save the PS mask for each cutout

pathOut = "./output/sehgal_maps/ir_sources/cutouts/"


#for iPatch in range(1,nPatches+1):


def generateMask(iPatch):
   print("Generating mask for cutout number "+str(iPatch))

   # read the cutout
   pathCutout = pathOut + "ir_"+sourceCatalog+"_sehgal_T_patch"+str(nPatches)+".txt"
   cutTMap = np.genfromtxt(pathCutout)

   # generate point source mask
   print("flux cut = "+str(fluxCutmJy)+" mJy")
   fluxCut = fluxCutmJy * 1.e-3 * 1.e-26   # convert from [mJy] to [Jy] to [W/m^2/Hz]
   fluxCut /= cmb.dBdT(nu, cmb.Tcmb)   # convert from [W/m^2/Hz] to [Kcmb*sr]
   fluxCut *= 1.e6   # convert from [Kcmb*sr] to [muKcmb*sr]
   print("ie flux cut = "+str(fluxCut)+" muKcmb*sr")

   # select patch radius around point sources
   maskPatchRadius = 3. * np.pi/(180.*60.)   # [arcmin] to [rad]
   
   # generate point source mask 
   cutTMapFourier = baseMap.fourier(cutTMap)
   psMask = baseMap.pointSourceMaskMatchedFilterIsotropic(cmb.ftotalTT, fluxCut, fprof=None, dataFourier=cutTMapFourier, maskPatchRadius=maskPatchRadius, test=False)    

   # save mask
   pathMask = pathOut + "ps_mask_ir_"+sourceCatalog+"_"+nuStr+"ghz_"+str(np.int(round(fluxCutmJy)))+"mJy_beam"+str(round(beamFwhm,1))+"_noise"+str(round(noiseT,2))+"_lknee"+str(np.int(lKnee))+"_aknee"+str(round(aKnee,1))+"_T_patch"+str(iPatch)+".txt"
   np.savetxt(pathMask, psMask)



nPatches = 41

# Cori Haswell as 32 cores per node, each accepting 2 threads
nProc = 41
with sharedmem.MapReduce(np=nProc) as pool:
   pool.map(generateMask, range(1, nPatches+1))

print("Done generating "+str(nPatches)+" masks")














