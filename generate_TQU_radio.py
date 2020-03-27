#!/usr/bin/env python
# coding: utf-8

# In[1]:


import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from basic_functions import *
from cmb import *
from flat_map import *


# # Read the radio galaxy catalog

# In[2]:


# path = "./input/sehgal_radio_catalog/radio_short.cat"
path = "./input/sehgal_radio_catalog/radio.cat"

#!!!! Just to make debugging quick, remove for final run!
# nObjMax = np.int(1.e7)

data = np.genfromtxt(path)#, max_rows=nObjMax)
ra = data[:,0]  # [deg]
dec = data[:,1]  # [deg]
z = data[:,2]
# f1_4_mJy = data[:,3]  # flux in mJy 
# f30_mJy = data[:,4]  # flux in mJy 
# f90_mJy = data[:,5]  # flux in mJy 
f148_mJy = data[:,6]  # flux in mJy 
# f219_mJy = data[:,7]  # flux in mJy 
# f277_mJy = data[:,8]  # flux in mJy 
# f350_mJy = data[:,9]  # flux in mJy 

print len(ra), "sources"


# Throw out the objects outside of the quadrant

# In[5]:


I = np.where((ra>=0.)*(ra<=90.)*(dec>=0.)*(dec<=90.))[0]
print "keeping", len(I), "objects out of", len(ra)
print "ie a fraction", 1.*len(I)/len(ra)


# In[6]:


ra = ra[I]
dec = dec[I]
z = z[I]
f148_mJy = f148_mJy[I]


# # Generate source count map (test)

# In[9]:


# Map geometry to match the Sehgal maps
nSide = 4096 #512#4096
nPix = hp.nside2npix(nSide)

# get pixel indices for all galaxies
IPix = hp.ang2pix(nSide, np.pi/2. - dec*np.pi/180., ra*np.pi/180., lonlat=False)

# print np.min(IPix), np.max(IPix)
# print nPix
# print IPix


# In[11]:


# Generate count map
bins = np.arange(nPix+1)-0.5
countMap, binEdges, binIndices = stats.binned_statistic(IPix, f148_mJy, statistic='count', bins=bins)

# IPix = 2*np.ones(nPix)
# countMap, binEdges = np.histogram(IPix, bins=nPix, range=(-0.5,nPix-0.5))


# # Generate T, Q, U maps

# In[15]:


# Generate T map
bins = np.arange(nPix+1)-0.5
tMap, binEdges, binIndices = stats.binned_statistic(IPix, f148_mJy, statistic='sum', bins=bins)  # flux map [mJy]

print "check that the map contains the flux from all the sources", np.sum(tMap), np.sum(f148_mJy)
print "ratio is", np.sum(tMap) / np.sum(f148_mJy)

tMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]


# In[16]:


# polarization fraction: 3% from Trombetti+18 (to be improved)
alpha = 0.03
# polarization angles:
theta = np.random.uniform(low=0., high=np.pi, size=len(ra))

# Generate Q and U maps
qMap, binEdges, binIndices = stats.binned_statistic(IPix, f148_mJy * alpha * np.cos(2.*theta), statistic='sum', bins=bins)  # flux map [mJy]
qMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]

uMap, binEdges, binIndices = stats.binned_statistic(IPix, f148_mJy * alpha * np.sin(2.*theta), statistic='sum', bins=bins)  # flux map [mJy]
uMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]


# # Compare T map to official Sehgal radio PS map

# In[17]:


# Read the official Sehgal map in temperature [Jy/sr]
sehgalTMap = hp.read_map("./input/sehgal_maps/148_rad_pts_healpix.fits")
# convert from [Jy/sr] to [mJy/sr]
sehgalTMap *= 1.e3


# In[18]:


sehgalTMap = hp.ud_grade(sehgalTMap, nSide)


# ## Convert all maps from [mJy/sr] to [muKcmb]

# In[21]:


cmb = CMB(beam=1., noise=1., nu1=148.e9, nu2=148.e9, lMin=30., lMaxT=3.e3, lMaxP=5.e3, fg=True, atm=False, name=None)

sehgalTMap *= 1.e-3 * 1.e-26  # convert from [mJy/sr] to surf bright per unit freq = [W/m^2/Hz/sr]
sehgalTMap /= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from surf bright per unit freq = [W/m^2/sr/Hz] to Kcmb
sehgalTMap *= 1.e6  # convert from Kcmb to muKcmb

tMap *= 1.e-3 * 1.e-26  # convert from [mJy/sr] to surf bright per unit freq = [W/m^2/Hz/sr]
tMap /= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from surf bright per unit freq = [W/m^2/sr/Hz] to Kcmb
tMap *= 1.e6  # convert from Kcmb to muKcmb

qMap *= 1.e-3 * 1.e-26  # convert from [mJy/sr] to surf bright per unit freq = [W/m^2/Hz/sr]
qMap /= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from surf bright per unit freq = [W/m^2/sr/Hz] to Kcmb
qMap *= 1.e6  # convert from Kcmb to muKcmb

uMap *= 1.e-3 * 1.e-26  # convert from [mJy/sr] to surf bright per unit freq = [W/m^2/Hz/sr]
uMap /= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from surf bright per unit freq = [W/m^2/sr/Hz] to Kcmb
uMap *= 1.e6  # convert from Kcmb to muKcmb


# The Lambda website says:
# dT = [Jy/sr] * T_CMB / 1.072480e9 in [T_CMB units]
# Check that it works:
print "My conversion agrees with the Lambda website recommendation:", 1.e-26 / cmb.dBdT(148.e9, cmb.Tcmb), cmb.Tcmb / 1.072480e9


# ## Read the kappa map

# In[24]:


kappaMap = hp.read_map("./input/sehgal_maps/healpix_4096_KappaeffLSStoCMBfullsky.fits")
kappaMap = hp.ud_grade(kappaMap, nSide)


# # Extract as many square cutouts as possible

# In[42]:


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


# In[ ]:


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
      if patchFits:
         nPatches += 1
         # plot the footprint
         xyz = hp.ang2vec(lonEdges, latEdges, lonlat=True)
         I = hp.query_polygon(nSide, xyz)
         hMap[I] += 1.
         
         # extract the cutouts
         pos = np.array([lonCenter, latCenter, 0.])
         # Official T map
         cutSehgalTMap = hp.visufunc.cartview(sehgalTMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.clf()
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
         cutKappaMap = hp.visufunc.cartview(kappaMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
         plt.clf()

         # generate point source mask
#          fluxCut = 5.1e-6  # for 2mJy
#         fluxCut = 12.75e-6  # for 5mJy
         fluxCut = 25.5e-6  # for 10mJy
         # convert to mJy to check
         fluxCutmJy = fluxCut * 1.e-6  # convert from [muKcmb*sr] to [Kcmb*sr]
         fluxCutmJy *= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from [Kcmb*sr] to flux per unit freq = [W/m^2/Hz]
         fluxCutmJy /= 1.e-26  # convert from flux per unit freq = [W/m^2/Hz] to [Jy]
         fluxCutmJy *= 1.e3  # convert from [Jy] to [mJy]
         print "ie", fluxCutmJy, "mJy"
         #
         # select patch radius around point sources
         maskPatchRadius = 3. * np.pi/(180.*60.)   # [arcmin] to [rad]
         #
         # generate point source mask 
         cutTMapFourier = baseMap.fourier(cutTMap)
         psMask = baseMap.pointSourceMaskMatchedFilterIsotropic(cmb.ftotalTT, fluxCut, fprof=None, dataFourier=cutTMapFourier, maskPatchRadius=maskPatchRadius, test=False)    
            
         # save the cutouts
         np.savetxt("./output/sehgal_maps/cutouts/ps_official_sehgal_T_patch"+str(nPatches)+".txt", cutSehgalTMap)
         np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_T_patch"+str(nPatches)+".txt", cutTMap)
         np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_Q_patch"+str(nPatches)+".txt", cutQMap)
         np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_U_patch"+str(nPatches)+".txt", cutUMap)
         np.savetxt("./output/sehgal_maps/cutouts/kappa_sehgal_patch"+str(nPatches)+".txt", cutKappaMap)
         np.savetxt("./output/sehgal_maps/cutouts/ps_mask_"+str(np.int(round(fluxCutmJy)))+"mJy_T_patch"+str(nPatches)+".txt", psMask)
         
print "Extracted "+str(nPatches)+" cutouts"




