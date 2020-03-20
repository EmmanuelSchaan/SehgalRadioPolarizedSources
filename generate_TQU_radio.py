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


# ## Check the location of the sources

# It is known that the Sehgal sims only have one octant. Weirdly, they are not quite exactly confined in the octant: some objects are slightly outside:

# In[3]:


print np.min(ra), np.max(ra)
print np.min(dec), np.max(dec)


# For an isotropic distribution of sources,
# the RA distribution should be flat
# and the dec distribution should go as $\cos(\text{dec})$

# In[4]:


#myHistogram(ra, nBins=71, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'RA', semilogx=False, semilogy=False, doGauss=False)
#myHistogram(dec, nBins=71, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'dec', semilogx=False, semilogy=False, doGauss=False)


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


# ## Check dn/dz of radio sources

# In[7]:


#myHistogram(z, nBins=71, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'$z$', semilogx=False, semilogy=False, doGauss=False)


# ## Check flux distribution

# The flux cut when masking will be around 5mJy.

# In[8]:


#myHistogram(f148_mJy, nBins=71, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'flux at 148GHz [mJy]', semilogx=True, semilogy=True, doGauss=False)


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


# The histogram of pixel indices should be flat, since the pixels have equal area and the sources are unclustered.

# In[10]:


#myHistogram(IPix, nBins=71, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'iPix', semilogx=False, semilogy=False, doGauss=False)


# In[11]:


# Generate count map
bins = np.arange(nPix+1)-0.5
countMap, binEdges, binIndices = stats.binned_statistic(IPix, f148_mJy, statistic='count', bins=bins)

# IPix = 2*np.ones(nPix)
# countMap, binEdges = np.histogram(IPix, bins=nPix, range=(-0.5,nPix-0.5))


# In[12]:


# print "bin edges", binEdges

print "check I have all objects", np.sum(countMap), len(ra)
print "check mean number per pixel", np.mean(countMap), 1. * len(ra) / nPix

print "check min / max number of objects per pixel", np.min(countMap), np.max(countMap)


# The pixel histogram should be a Poisson distribution

# In[13]:


#myHistogram(1.e-10+countMap, nBins=101, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'counts', semilogx=False, semilogy=True, doGauss=False)


# Check that the correct quadrant in the sky is filled: $\text{RA}, \text{dec}\in[0, \pi/2]$.

# In[14]:


#hp.mollview(np.log10(1.e-1+countMap))


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


# In[ ]:


sehgalTMap = hp.ud_grade(sehgalTMap, nSide)


# In[ ]:


#myHistogram(1.e-10 + tMap, nBins=501, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'T', semilogx=True, semilogy=True, doGauss=False)
#myHistogram(1.e-10 + sehgalTMap, nBins=501, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'T', semilogx=True, semilogy=True, doGauss=False)

#myHistogram(1.e-10 + np.abs(tMap - sehgalTMap), nBins=501, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'T', semilogx=True, semilogy=True, doGauss=False)


# In[ ]:


#hp.mollview(np.log10(1.e-1+sehgalTMap), min=-1, max=10)
#hp.mollview(np.log10(1.e-1+tMap), min=-1, max=10)
#hp.mollview(np.log10(1.e-10 + np.abs(sehgalTMap - tMap)), min=np.log10(1.e-5), max=np.log10(1.e10))


# ## Convert all maps from [mJy/sr] to [muKcmb]

# In[ ]:


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


# In[ ]:


#myHistogram(1.e-4+tMap, nBins=71, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'$T$ [$\mu$K$_\text{CMB}$]', semilogx=True, semilogy=True, doGauss=False)
#myHistogram(1.e-11+np.abs(qMap), nBins=71, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'$|Q|$ [$\mu$K$_\text{CMB}$]', semilogx=True, semilogy=True, doGauss=False)
#myHistogram(1.e-11+np.abs(uMap), nBins=71, lim=None, S2Theory=[], path=None, plot=True, nameLatex=r'$|U|$ [$\mu$K$_\text{CMB}$]', semilogx=True, semilogy=True, doGauss=False)


# In[ ]:


#hp.mollview(np.log10(1.e-4+tMap))
#hp.mollview(np.log10(1.e-11+qMap))
#hp.mollview(np.log10(1.e-11+uMap))


# ## Read the kappa map

# In[ ]:


kappaMap = hp.read_map("./input/sehgal_maps/healpix_4096_KappaeffLSStoCMBfullsky.fits")
kappaMap = hp.ud_grade(kappaMap, nSide)


# # Extract cutouts T, Q, U, official Sehgal T (test) and kappa

# In[ ]:


# position (lon, lat, psi) of the cutout center
lon = 45. # [deg]
lat = 45. # [deg]
pos = np.array([lon, lat, 0.]) 

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

print xSize, "pixel on the side,", dLon, "degrees on the side"
print "resolution is", dLon/xSize*60., "arcmin"


# In[ ]:


# Official T map
# cutTMap = hp.visufunc.cartview(tMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
cutSehgalTMap = hp.visufunc.cartview(sehgalTMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
# T map
cutTMap = hp.visufunc.cartview(tMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
# Q map
cutQMap = hp.visufunc.cartview(qMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
# U map
cutUMap = hp.visufunc.cartview(uMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')
# kappa map
cutKappaMap = hp.visufunc.cartview(kappaMap, rot=pos, lonra=lonRange, latra=latRange, xsize=xSize, ysize=ySize, return_projected_map=True, norm='hist')


# ## Test: the cutout from my T map and the official T map should be the same.

# They have the same power spectrum:

# In[ ]:


cutTMapFourier = baseMap.fourier(cutTMap)
cutSehgalTMapFourier = baseMap.fourier(cutSehgalTMap)

lCen, Cl, sCl = baseMap.powerSpectrum(cutTMapFourier, plot=False)
lCen, ClSehgal, sCl = baseMap.powerSpectrum(cutSehgalTMapFourier, plot=False)
# lCen, Cl, sCl = baseMap.powerSpectrum(cutTMapFourier - cutSehgalTMapFourier, plot=True)

plt.semilogx(lCen, Cl/ClSehgal)
#plt.show()


# The map difference has negligible power:

# In[ ]:


lCen, Cl, sCl = baseMap.powerSpectrum(cutTMapFourier, plot=False)
lCen, ClRes, sCl = baseMap.powerSpectrum(cutTMapFourier - cutSehgalTMapFourier, plot=False)

plt.loglog(lCen, ClRes / Cl)
#plt.show()


# ## Test: mask the point sources above 15mJy, and compare power to Dunkley+13

# In[ ]:


# use the flux cut from Dunkley+13, to match power spectra
fluxCut = 0.015  # in [Jy]

# convert to muKcmb
fluxCut *= 1.e-26  # convert from [Jy] to flux per unit freq = [W/m^2/Hz]
fluxCut /= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from flux per unit freq = [W/m^2/Hz] to [Kcmb*sr]
fluxCut *= 1.e6  # convert from [Kcmb*sr] to [muKcmb*sr]


# In[ ]:


maskPatchRadius = 3. * np.pi/(180.*60.)   # [arcmin] to [rad]
psMask = baseMap.pointSourceMaskMatchedFilterIsotropic(cmb.ftotalTT, fluxCut, fprof=None, dataFourier=cutTMapFourier, maskPatchRadius=maskPatchRadius, test=False)


# In[ ]:


#baseMap.plot(psMask)


# Power spectrum before masking: much higher than Dunkley+13

# In[ ]:


lCen, Cl, sCl = baseMap.powerSpectrum(cutTMapFourier, theory=[cmb.fradioPoisson], plot=True)


# Power spectrum after masking: close enough to Dunkley+13

# In[ ]:


maskedTMapFourier = baseMap.fourier(cutTMap*psMask)
lCen, Cl, sCl = baseMap.powerSpectrum(maskedTMapFourier, theory=[cmb.fradioPoisson], plot=True)


# # Use point source mask relevant for CMB S4

# In[ ]:


# experimental specs for CMB S4
cmb = CMB(beam=1., noise=1., nu1=148.e9, nu2=148.e9, lMin=30., lMaxT=1.e5, lMaxP=1.e5, fg=True, atm=False, name=None)

# find the 5 sigma point source detection threshold [muK*sr]
fluxCut = 5. * cmb.fsigmaMatchedFilter(fprofile=None, ftotalTT=cmb.ftotalTT, lMin=30., lMax=1.e4)

print "the 5 sigma flux cut is", fluxCut, "muK*sr"


# Round it off to get a 2mJy point source cut

# In[ ]:


# rounded off flux cut to get 2mJy
fluxCut = 5.1e-6
# convert to mJy to check
fluxCutmJy = fluxCut * 1.e-6  # convert from [muKcmb*sr] to [Kcmb*sr]
fluxCutmJy *= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from [Kcmb*sr] to flux per unit freq = [W/m^2/Hz]
fluxCutmJy /= 1.e-26  # convert from flux per unit freq = [W/m^2/Hz] to [Jy]
fluxCutmJy *= 1.e3  # convert from [Jy] to [mJy]
print "ie", fluxCutmJy, "mJy"

# select patch radius around point sources
maskPatchRadius = 3. * np.pi/(180.*60.)   # [arcmin] to [rad]

# generate point source mask 
psMask = baseMap.pointSourceMaskMatchedFilterIsotropic(cmb.ftotalTT, fluxCut, fprof=None, dataFourier=cutTMapFourier, maskPatchRadius=maskPatchRadius, test=False)


# In[ ]:


#baseMap.plot(psMask)


# # Save all cutout maps and mask

# In[ ]:


np.savetxt("./output/sehgal_maps/cutouts/x_rad.txt", baseMap.x)
np.savetxt("./output/sehgal_maps/cutouts/y_rad.txt", baseMap.y)
np.savetxt("./output/sehgal_maps/cutouts/ps_official_sehgal_T.txt", cutSehgalTMap)
np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_T.txt", cutTMap)
np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_Q.txt", cutQMap)
np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_U.txt", cutUMap)
np.savetxt("./output/sehgal_maps/cutouts/kappa_sehgal.txt", cutKappaMap)
np.savetxt("./output/sehgal_maps/cutouts/ps_mask_2mJy_T.txt", psMask)


# # Redo it on a different patch

# In[ ]:


# position (lon, lat, psi) of the cutout center
lon = 25. # [deg]
lat = 25. # [deg]
pos = np.array([lon, lat, 0.]) 

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

print xSize, "pixel on the side,", dLon, "degrees on the side"
print "resolution is", dLon/xSize*60., "arcmin"


# In[ ]:


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


# In[ ]:


# rounded off flux cut to get 2mJy
fluxCut = 5.1e-6
# convert to mJy to check
fluxCutmJy = fluxCut * 1.e-6  # convert from [muKcmb*sr] to [Kcmb*sr]
fluxCutmJy *= cmb.dBdT(148.e9, cmb.Tcmb)  # convert from [Kcmb*sr] to flux per unit freq = [W/m^2/Hz]
fluxCutmJy /= 1.e-26  # convert from flux per unit freq = [W/m^2/Hz] to [Jy]
fluxCutmJy *= 1.e3  # convert from [Jy] to [mJy]
print "ie", fluxCutmJy, "mJy"

# select patch radius around point sources
maskPatchRadius = 3. * np.pi/(180.*60.)   # [arcmin] to [rad]

# generate point source mask 
cutTMapFourier = baseMap.fourier(cutTMap)
psMask = baseMap.pointSourceMaskMatchedFilterIsotropic(cmb.ftotalTT, fluxCut, fprof=None, dataFourier=cutTMapFourier, maskPatchRadius=maskPatchRadius, test=False)


# In[ ]:


np.savetxt("./output/sehgal_maps/cutouts/ps_official_sehgal_T_2.txt", cutSehgalTMap)
np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_T_2.txt", cutTMap)
np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_Q_2.txt", cutQMap)
np.savetxt("./output/sehgal_maps/cutouts/ps_sehgal_U_2.txt", cutUMap)
np.savetxt("./output/sehgal_maps/cutouts/kappa_sehgal_2.txt", cutKappaMap)
np.savetxt("./output/sehgal_maps/cutouts/ps_mask_2mJy_T_2.txt", psMask)


# # Extract as many square cutouts as possible

# In[ ]:


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
         # rounded off flux cut to get 2mJy
         fluxCut = 5.1e-6
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
         np.savetxt("./output/sehgal_maps/cutouts/ps_mask_2mJy_T_patch"+str(nPatches)+".txt", psMask)
         
print "Extracted "+str(nPatches)+" cutouts"

#hp.mollview(hMap)
#plt.show()


# In[ ]:




