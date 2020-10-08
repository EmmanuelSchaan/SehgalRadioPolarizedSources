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


nu = 90.e9  # [Hz]
nuStr = '90'
fluxCutmJy = 2.  # [mJy]

#nu = 148.e9  # [Hz]
#nuStr = '148'
#fluxCutmJy = 2.  # [mJy]


#####################################################################
## # Read the radio galaxy catalog
#
## path = "./input/sehgal_radio_catalog/radio_short.cat"
#path = "./input/sehgal_radio_catalog/radio.cat"
#
## Just to make debugging quick, remove for final run
## nObjMax = np.int(1.e7)
#
#data = np.genfromtxt(path)#, max_rows=nObjMax)
#ra = data[:,0]  # [deg]
#dec = data[:,1]  # [deg]
#z = data[:,2]
#
#fluxes_mJy = {}
#fluxes_mJy['1_4'] = data[:,3]  # flux in mJy 
#fluxes_mJy['30'] = data[:,4]  # flux in mJy 
#fluxes_mJy['90'] = data[:,5]  # flux in mJy 
#fluxes_mJy['148'] = data[:,6]  # flux in mJy 
#fluxes_mJy['219'] = data[:,7]  # flux in mJy 
#fluxes_mJy['277'] = data[:,8]  # flux in mJy 
#fluxes_mJy['350'] = data[:,9]  # flux in mJy 
#
## select the requested frequency
#flux_mJy = fluxes_mJy[nuStr]
#
#print len(ra), "sources"
#
#
## Throw out the objects outside of the quadrant
#I = np.where((ra>=0.)*(ra<=90.)*(dec>=0.)*(dec<=90.))[0]
#print "keeping", len(I), "objects out of", len(ra)
#print "ie a fraction", 1.*len(I)/len(ra)
#
#ra = ra[I]
#dec = dec[I]
#z = z[I]
#flux_mJy = flux_mJy[I]
#
#
#####################################################################
## # Generate source count map (for testing purposes)
#
## Map geometry to match the Sehgal maps
#nSide = 4096 #512#4096
#nPix = hp.nside2npix(nSide)
#
## get pixel indices for all galaxies
#IPix = hp.ang2pix(nSide, np.pi/2. - dec*np.pi/180., ra*np.pi/180., lonlat=False)
#
## Generate count map
#bins = np.arange(nPix+1)-0.5
#countMap, binEdges, binIndices = stats.binned_statistic(IPix, flux_mJy, statistic='count', bins=bins)
#
#
#####################################################################
## # Generate T, Q, U maps
#
## Generate T map
#bins = np.arange(nPix+1)-0.5
#tMap, binEdges, binIndices = stats.binned_statistic(IPix, flux_mJy, statistic='sum', bins=bins)  # flux map [mJy]
#
#print "check that the map contains the flux from all the sources", np.sum(tMap), np.sum(flux_mJy)
#print "ratio is", np.sum(tMap) / np.sum(flux_mJy)
#
#tMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]
#
#
## polarization fraction: 3% from Trombetti+18 (to be improved)
#alpha = 0.03
## polarization angles:
#theta = np.random.uniform(low=0., high=np.pi, size=len(ra))
#
## Generate Q and U maps
#qMap, binEdges, binIndices = stats.binned_statistic(IPix, flux_mJy * alpha * np.cos(2.*theta), statistic='sum', bins=bins)  # flux map [mJy]
#qMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]
#
#uMap, binEdges, binIndices = stats.binned_statistic(IPix, flux_mJy * alpha * np.sin(2.*theta), statistic='sum', bins=bins)  # flux map [mJy]
#uMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]
#
#
#####################################################################
## # Compare T map to official Sehgal radio PS map
#
## Read the official Sehgal map in temperature [Jy/sr]
##sehgalTMap = hp.read_map("./input/sehgal_maps/148_rad_pts_healpix.fits")
#sehgalTMap = hp.read_map("./input/sehgal_maps/"+nuStr.zfill(3)+"_rad_pts_healpix.fits")
## convert from [Jy/sr] to [mJy/sr]
#sehgalTMap *= 1.e3
#sehgalTMap = hp.ud_grade(sehgalTMap, nSide)
#
#
## ## Convert all maps from [mJy/sr] to [muKcmb]
#cmb = CMB(beam=1., noise=1., nu1=nu, nu2=nu, lMin=30., lMaxT=3.e3, lMaxP=5.e3, fg=True, atm=False, name=None)
#
#conversionFactor = 1.e-3 * 1.e-26  # convert from [mJy/sr] to surf bright per unit freq = [W/m^2/Hz/sr]
#conversionFactor /= cmb.dBdT(nu, cmb.Tcmb)  # convert from surf bright per unit freq = [W/m^2/sr/Hz] to Kcmb
#conversionFactor *= 1.e6  # convert from Kcmb to muKcmb
#
#sehgalTMap *= conversionFactor
#tMap *= conversionFactor
#qMap *= conversionFactor
#uMap *= conversionFactor
#
#
## The Lambda website says:
## dT = [Jy/sr] * T_CMB / conversionOfficial in [T_CMB units]
#conversionOfficial = {}
#conversionOfficial['30'] = 7.364967e7
#conversionOfficial['90'] = 5.526540e8
#conversionOfficial['148'] = 1.072480e9
#conversionOfficial['219'] = 1.318837e9
#conversionOfficial['277'] = 1.182877e9
#conversionOfficial['350'] = 8.247628e8
## Check that it works:
#print "My conversion agrees with the Lambda website recommendation at "+nuStr+" GHz:", 1.e-26 / cmb.dBdT(nu, cmb.Tcmb), cmb.Tcmb / conversionOfficial[nuStr]
#
#
#####################################################################
#print "Saving healpix maps"
#
##pathOut = "/global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/output/sehgal_maps/radio_sources/"
#pathOut = "./output/sehgal_maps/radio_sources/"
#hp.write_map(pathOut + "t_radio_sehgal_"+nuStr+"ghz_muk.fits", tMap, overwrite=True)
#hp.write_map(pathOut + "q_radio_sehgal_"+nuStr+"ghz_muk.fits", qMap, overwrite=True)
#hp.write_map(pathOut + "u_radio_sehgal_"+nuStr+"ghz_muk.fits", uMap, overwrite=True)


#####################################################################
print "Reading healpix maps"

sehgalTMap = hp.read_map("./input/sehgal_maps/"+nuStr.zfill(3)+"_rad_pts_healpix.fits")

#pathOut = "/global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/output/sehgal_maps/radio_sources/"
pathOut = "./output/sehgal_maps/radio_sources/"
tMap = hp.read_map(pathOut + "t_radio_sehgal_"+nuStr+"ghz_muk.fits")
qMap = hp.read_map(pathOut + "q_radio_sehgal_"+nuStr+"ghz_muk.fits")
uMap = hp.read_map(pathOut + "u_radio_sehgal_"+nuStr+"ghz_muk.fits")


nSide = hp.get_nside(tMap)

cmb = CMB(beam=1., noise=1., nu1=nu, nu2=nu, lMin=30., lMaxT=3.e3, lMaxP=5.e3, fg=True, atm=False, name=None)

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
#          fluxCut = 5.1e-6  # for 2mJy
#         fluxCut = 12.75e-6  # for 5mJy
#         fluxCut = 25.5e-6  # for 10mJy
#         # convert to mJy to check
#         fluxCutmJy = fluxCut * 1.e-6  # convert from [muKcmb*sr] to [Kcmb*sr]
#         fluxCutmJy *= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from [Kcmb*sr] to flux per unit freq = [W/m^2/Hz]
#         fluxCutmJy /= 1.e-26  # convert from flux per unit freq = [W/m^2/Hz] to [Jy]
#         fluxCutmJy *= 1.e3  # convert from [Jy] to [mJy]
#         print "ie", fluxCutmJy, "mJy"
         print("flux cut = "+str(fluxCutmJy)+" mJy")
         fluxCut = fluxCutmJy * 1.e-3 * 1.e-26   # convert from [mJy] to [Jy] to [W/m^2/Hz]
         fluxCut /= cmb.dBdT(nu, cmb.Tcmb)   # convert from [W/m^2/Hz] to [Kcmb*sr]
         fluxCut *= 1.e6   # convert from [Kcmb*sr] to [muKcmb*sr]
         print("ie flux cut = "+str(fluxCut)+" muKcmb*sr")



         #
         # select patch radius around point sources
         maskPatchRadius = 3. * np.pi/(180.*60.)   # [arcmin] to [rad]
         #
         # generate point source mask 
         cutTMapFourier = baseMap.fourier(cutTMap)
         psMask = baseMap.pointSourceMaskMatchedFilterIsotropic(cmb.ftotalTT, fluxCut, fprof=None, dataFourier=cutTMapFourier, maskPatchRadius=maskPatchRadius, test=False)    
            
#         # save the cutouts
#         np.savetxt("./output/sehgal_maps/cutouts/ps_official_sehgal_T_patch"+str(nPatches)+".txt", cutSehgalTMap)
#         np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_T_patch"+str(nPatches)+".txt", cutTMap)
#         np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_Q_patch"+str(nPatches)+".txt", cutQMap)
#         np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_U_patch"+str(nPatches)+".txt", cutUMap)
#         np.savetxt("./output/sehgal_maps/cutouts/kappa_sehgal_patch"+str(nPatches)+".txt", cutKappaMap)
#         np.savetxt("./output/sehgal_maps/cutouts/ps_mask_"+str(np.int(round(fluxCutmJy)))+"mJy_T_patch"+str(nPatches)+".txt", psMask)

         # save the cutouts
         #pathOut = "/global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/output/sehgal_maps/radio_sources/cutouts/"
         pathOut = "./output/sehgal_maps/radio_sources/cutouts/"
         np.savetxt(pathOut + "ps_official_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt", cutSehgalTMap)
         np.savetxt(pathOut + "ps_sehgal_"+nuStr+"ghz_T_patch"+str(nPatches)+".txt", cutTMap)
         np.savetxt(pathOut + "ps_sehgal_"+nuStr+"ghz_Q_patch"+str(nPatches)+".txt", cutQMap)
         np.savetxt(pathOut + "ps_sehgal_"+nuStr+"ghz_U_patch"+str(nPatches)+".txt", cutUMap)
         np.savetxt(pathOut + "kappa_sehgal_patch"+str(nPatches)+".txt", cutKappaMap)
         np.savetxt(pathOut + "ps_mask_"+nuStr+"ghz_"+str(np.int(round(fluxCutmJy)))+"mJy_T_patch"+str(nPatches)+".txt", psMask)

         
print "Extracted "+str(nPatches)+" cutouts"




