import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from basic_functions import *
from cmb import *
from flat_map import *

import sys


# needed on lawrencium
plt.switch_backend('Agg')

####################################################################

# # Read one of the galaxy catalog
# Choose input file
sourceCatalog = sys.argv[1]

print("Analyzing catalog "+sourceCatalog)

####################################################################

pathIn = "./output/sehgal_maps/ir_sources/"
pathOut = "./output/sehgal_maps/ir_sources/cutouts/"

# read healpix maps
tMap = hp.read_map(pathIn + "t_ir_"+sourceCatalog+"_sehgal_148ghz_muk.fits")
qMap = hp.read_map(pathIn + "q_ir_"+sourceCatalog+"_sehgal_148ghz_muk.fits")
uMap = hp.read_map(pathIn + "u_ir_"+sourceCatalog+"_sehgal_148ghz_muk.fits")
print("Done reading the T, Q, U maps")

#cmb = CMB(beam=1., noise=1., nu1=148.e9, nu2=148.e9, lMin=30., lMaxT=3.e3, lMaxP=5.e3, fg=True, atm=False, name=None)

print("Extract as many square cutouts as possible")

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

# create an empty map, to check footprints
nSide = 256
hMap = np.zeros(hp.nside2npix(nSide))


# offsets in lon and lat to start the cutouts
latStart = 0. #1.
lonStart = 0. #1.
# space between cutouts, to avoid overlap
space = 0.5 # [deg]


# latitudes of centers of cutouts
iPatch = 0
LatCenter = np.arange(latStart + dLat/2., 90., dLat + space)
for latCenter in LatCenter:
   latUpper = latCenter + dLat/2.
   latLower = latCenter - dLat/2.
   
   dLonCenter = dLon / np.cos(latCenter * np.pi/180.)
   dLonUpper = dLon / np.cos(latUpper * np.pi/180.)
   dLonLower = dLon / np.cos(latLower * np.pi/180.)
   LonCenter = np.arange(lonStart + dLonUpper/2., 90., dLonUpper + space)
   
   for lonCenter in LonCenter:
      
      # polygon edges
      latEdges = np.array([latUpper, latLower, latLower, latUpper])
      lonEdges = np.array([lonCenter-dLonUpper/2., lonCenter-dLonLower/2., lonCenter+dLonLower/2., lonCenter+dLonUpper/2.])
      
      # check that the cutout is entirely within the quadrant
      patchFits = np.sum(latEdges<0.) + np.sum(latEdges>90.) + np.sum(lonEdges<0.) + np.sum(lonEdges>90.)
      patchFits = patchFits==0.
      
      # if it fits
      if patchFits:
         iPatch += 1
         print("Extracting cutout number "+str(iPatch))

         # plot the footprint
         xyz = hp.ang2vec(lonEdges, latEdges, lonlat=True)
         I = hp.query_polygon(nSide, xyz)
         hMap[I] += 1.
         
         # extract the cutouts
         pos = np.array([lonCenter, latCenter, 0.])
         # Official T map
#          cutSehgalTMap = hp.visufunc.cartview(sehgalTMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
#          plt.clf()
         # T map
         cutTMap = hp.visufunc.cartview(tMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.clf()
         # Q map
         cutQMap = hp.visufunc.cartview(qMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.clf()
         # U map
         cutUMap = hp.visufunc.cartview(uMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.clf()
         # kappa map
#          cutKappaMap = hp.visufunc.cartview(kappaMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
#          plt.clf()

#         # generate point source mask
#         # fluxCut = 5.1e-6  # for 2mJy
#         fluxCut = 12.75e-6  # for 5mJy
#         # convert to mJy to check
#         fluxCutmJy = fluxCut * 1.e-6  # convert from [muKcmb*sr] to [Kcmb*sr]
#         fluxCutmJy *= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from [Kcmb*sr] to flux per unit freq = [W/m^2/Hz]
#         fluxCutmJy /= 1.e-26  # convert from flux per unit freq = [W/m^2/Hz] to [Jy]
#         fluxCutmJy *= 1.e3  # convert from [Jy] to [mJy]
#         print("ie", fluxCutmJy, "mJy")
#         #
#         # select patch radius around point sources
#         maskPatchRadius = 3. * np.pi/(180.*60.)   # [arcmin] to [rad]
#         #
#         # generate point source mask 
#         cutTMapFourier = baseMap.fourier(cutTMap)
#         psMask = baseMap.pointSourceMaskMatchedFilterIsotropic(cmb.ftotalTT, fluxCut, fprof=None, dataFourier=cutTMapFourier, maskPatchRadius=maskPatchRadius, test=False)    
            
         # save the cutouts
#          np.savetxt("./output/sehgal_maps/cutouts/ir_official_sehgal_T_patch"+str(iPatch)+".txt", cutSehgalTMap)
         np.savetxt(pathOut + "ir_"+sourceCatalog+"_sehgal_T_patch"+str(iPatch)+".txt", cutTMap)
         np.savetxt(pathOut + "ir_"+sourceCatalog+"_sehgal_Q_patch"+str(iPatch)+".txt", cutQMap)
         np.savetxt(pathOut + "ir_"+sourceCatalog+"_sehgal_U_patch"+str(iPatch)+".txt", cutUMap)
#          np.savetxt("./output/sehgal_maps/cutouts/kappa_sehgal_patch"+str(iPatch)+".txt", cutKappaMap)
#         np.savetxt(pathOut + "ir_"+sourceCatalog+"_mask_"+str(np.int(round(fluxCutmJy)))+"mJy_T_patch"+str(iPatch)+".txt", psMask)
         
print("All finished!")
print("Extracted "+str(iPatch)+" cutouts")






