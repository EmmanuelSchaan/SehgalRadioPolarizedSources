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




#####################################################################
print "Reading healpix maps"

sehgalTMap = hp.read_map("./input/sehgal_maps/"+nuStr.zfill(3)+"_rad_pts_healpix.fits")

#pathOut = "/global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/output/sehgal_maps/radio_sources/"
pathOut = "./output/sehgal_maps/radio_sources/"
tMap = hp.read_map(pathOut + "t_radio_sehgal_"+nuStr+"ghz_muk.fits")
qMap = hp.read_map(pathOut + "q_radio_sehgal_"+nuStr+"ghz_muk.fits")
uMap = hp.read_map(pathOut + "u_radio_sehgal_"+nuStr+"ghz_muk.fits")


nSide = hp.get_nside(tMap)

cmb = CMB(beam=beamFwhm, noise=noiseT, nu1=nu, nu2=nu, lMin=30., lMaxT=3.e3, lMaxP=5.e3, fg=True, atm=True, atmProp=[lKnee, aKnee, lKnee, aKnee], name=None)

####################################################################
# ## Read the kappa map

kappaMap = hp.read_map("./input/sehgal_maps/healpix_4096_KappaeffLSStoCMBfullsky.fits")
kappaMap = hp.ud_grade(kappaMap, nSide)


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

# create an empty map, to check footprints
nSide = 256
hMap = np.zeros(hp.nside2npix(nSide))


# offsets in lon and lat to start the cutouts
latStart = 0. #1.
lonStart = 0. #1.
# space between cutouts, to avoid overlap
space = 0.5 # [deg]




# latitudes of centers of cutouts
nPatches = 0
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
      #if patchFits and nPatches>=17:
      if patchFits: 
         nPatches += 1

         #if nPatches>=17:
         #if nPatches>=30:

         # plot the footprint
         xyz = hp.ang2vec(lonEdges, latEdges, lonlat=True)
         I = hp.query_polygon(nSide, xyz)
         hMap[I] += 1.
         
         # extract the cutouts
         pos = np.array([lonCenter, latCenter, 0.])
         # Official T map
         cutSehgalTMap = hp.visufunc.cartview(sehgalTMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.close()
         # T map
         cutTMap = hp.visufunc.cartview(tMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.close()
         # Q map
         cutQMap = hp.visufunc.cartview(qMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.close()
         # U map
         cutUMap = hp.visufunc.cartview(uMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.close()
         # kappa map
         cutKappaMap = hp.visufunc.cartview(kappaMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.close()

         # generate point source mask
         print("flux cut = "+str(fluxCutmJy)+" mJy")
         fluxCut = fluxCutmJy * 1.e-3 * 1.e-26   # convert from [mJy] to [Jy] to [W/m^2/Hz]
         fluxCut /= cmb.dBdT(nu, cmb.Tcmb)   # convert from [W/m^2/Hz] to [Kcmb*sr]
         fluxCut *= 1.e6   # convert from [Kcmb*sr] to [muKcmb*sr]
         print("ie flux cut = "+str(fluxCut)+" muKcmb*sr")


         # select patch radius around point sources
         maskPatchRadius = 3. * np.pi/(180.*60.)   # [arcmin] to [rad]
         #
         # generate point source mask 
         cutTMapFourier = baseMap.fourier(cutTMap)
         psMask = baseMap.pointSourceMaskMatchedFilterIsotropic(cmb.ftotalTT, fluxCut, fprof=None, dataFourier=cutTMapFourier, maskPatchRadius=maskPatchRadius, test=False)    
            
         # save the cutouts
         #pathOut = "/global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/output/sehgal_maps/radio_sources/cutouts/"
         pathOut = "./output/sehgal_maps/radio_sources/cutouts/"
#         np.savetxt(pathOut + "ps_official_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt", cutSehgalTMap)
#         np.savetxt(pathOut + "ps_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt", cutTMap)
#         np.savetxt(pathOut + "ps_sehgal_"+nuStr+"ghz_Q_patch"+str(nPatches)+".txt", cutQMap)
#         np.savetxt(pathOut + "ps_sehgal_"+nuStr+"ghz_U_patch"+str(nPatches)+".txt", cutUMap)
#         np.savetxt(pathOut + "kappa_sehgal_patch"+str(nPatches)+".txt", cutKappaMap)
         np.savetxt(pathOut + "ps_mask_"+nuStr+"ghz_"+str(np.int(round(fluxCutmJy)))+"mJy_beam"+str(round(beamFwhm,1))+"_noise"+str(round(noiseT,2))+"_lknee"+str(np.int(lKnee))+"_aknee"+str(round(aKnee,1))+"_T_patch"+str(nPatches)+".txt", psMask)

         
print "Extracted "+str(nPatches)+" cutouts"




