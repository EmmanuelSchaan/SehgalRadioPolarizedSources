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

#nu = 148.e9  # [Hz]
#nuStr = '148'


####################################################################
# # Read the radio galaxy catalog

# path = "./input/sehgal_radio_catalog/radio_short.cat"
path = "./input/sehgal_radio_catalog/radio.cat"

# Just to make debugging quick, remove for final run
# nObjMax = np.int(1.e7)

data = np.genfromtxt(path)#, max_rows=nObjMax)
ra = data[:,0]  # [deg]
dec = data[:,1]  # [deg]
z = data[:,2]

fluxes_mJy = {}
fluxes_mJy['1_4'] = data[:,3]  # flux in mJy 
fluxes_mJy['30'] = data[:,4]  # flux in mJy 
fluxes_mJy['90'] = data[:,5]  # flux in mJy 
fluxes_mJy['148'] = data[:,6]  # flux in mJy 
fluxes_mJy['219'] = data[:,7]  # flux in mJy 
fluxes_mJy['277'] = data[:,8]  # flux in mJy 
fluxes_mJy['350'] = data[:,9]  # flux in mJy 

# select the requested frequency
flux_mJy = fluxes_mJy[nuStr]

print len(ra), "sources"


# Throw out the objects outside of the quadrant
I = np.where((ra>=0.)*(ra<=90.)*(dec>=0.)*(dec<=90.))[0]
print "keeping", len(I), "objects out of", len(ra)
print "ie a fraction", 1.*len(I)/len(ra)

ra = ra[I]
dec = dec[I]
z = z[I]
flux_mJy = flux_mJy[I]


####################################################################
# # Generate source count map (for testing purposes)

# Map geometry to match the Sehgal maps
nSide = 4096 #512#4096
nPix = hp.nside2npix(nSide)

# get pixel indices for all galaxies
IPix = hp.ang2pix(nSide, np.pi/2. - dec*np.pi/180., ra*np.pi/180., lonlat=False)

# Generate count map
bins = np.arange(nPix+1)-0.5
countMap, binEdges, binIndices = stats.binned_statistic(IPix, flux_mJy, statistic='count', bins=bins)


####################################################################
# # Generate T, Q, U maps

# Generate T map
bins = np.arange(nPix+1)-0.5
tMap, binEdges, binIndices = stats.binned_statistic(IPix, flux_mJy, statistic='sum', bins=bins)  # flux map [mJy]

print "check that the map contains the flux from all the sources", np.sum(tMap), np.sum(flux_mJy)
print "ratio is", np.sum(tMap) / np.sum(flux_mJy)

tMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]


# polarization fraction: 3% from Trombetti+18 (to be improved)
alpha = 0.03
# polarization angles:
theta = np.random.uniform(low=0., high=np.pi, size=len(ra))

# Generate Q and U maps
qMap, binEdges, binIndices = stats.binned_statistic(IPix, flux_mJy * alpha * np.cos(2.*theta), statistic='sum', bins=bins)  # flux map [mJy]
qMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]

uMap, binEdges, binIndices = stats.binned_statistic(IPix, flux_mJy * alpha * np.sin(2.*theta), statistic='sum', bins=bins)  # flux map [mJy]
uMap /= hp.nside2pixarea(nSide)  # surf bright map [mJy/sr]


####################################################################
# # Compare T map to official Sehgal radio PS map

# Read the official Sehgal map in temperature [Jy/sr]
#sehgalTMap = hp.read_map("./input/sehgal_maps/148_rad_pts_healpix.fits")
sehgalTMap = hp.read_map("./input/sehgal_maps/"+nuStr.zfill(3)+"_rad_pts_healpix.fits")
# convert from [Jy/sr] to [mJy/sr]
sehgalTMap *= 1.e3
sehgalTMap = hp.ud_grade(sehgalTMap, nSide)


# ## Convert all maps from [mJy/sr] to [muKcmb]
cmb = CMB(beam=1., noise=1., nu1=nu, nu2=nu, lMin=30., lMaxT=3.e3, lMaxP=5.e3, fg=True, atm=False, name=None)

conversionFactor = 1.e-3 * 1.e-26  # convert from [mJy/sr] to surf bright per unit freq = [W/m^2/Hz/sr]
conversionFactor /= cmb.dBdT(nu, cmb.Tcmb)  # convert from surf bright per unit freq = [W/m^2/sr/Hz] to Kcmb
conversionFactor *= 1.e6  # convert from Kcmb to muKcmb

sehgalTMap *= conversionFactor
tMap *= conversionFactor
qMap *= conversionFactor
uMap *= conversionFactor


# The Lambda website says:
# dT = [Jy/sr] * T_CMB / conversionOfficial in [T_CMB units]
conversionOfficial = {}
conversionOfficial['30'] = 7.364967e7
conversionOfficial['90'] = 5.526540e8
conversionOfficial['148'] = 1.072480e9
conversionOfficial['219'] = 1.318837e9
conversionOfficial['277'] = 1.182877e9
conversionOfficial['350'] = 8.247628e8
# Check that it works:
print "My conversion agrees with the Lambda website recommendation at "+nuStr+" GHz:", 1.e-26 / cmb.dBdT(nu, cmb.Tcmb), cmb.Tcmb / conversionOfficial[nuStr]


####################################################################
print "Saving healpix maps"

#pathOut = "/global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/output/sehgal_maps/radio_sources/"
pathOut = "./output/sehgal_maps/radio_sources/"
hp.write_map(pathOut + "t_radio_sehgal_"+nuStr+"ghz_muk.fits", tMap, overwrite=True)
hp.write_map(pathOut + "q_radio_sehgal_"+nuStr+"ghz_muk.fits", qMap, overwrite=True)
hp.write_map(pathOut + "u_radio_sehgal_"+nuStr+"ghz_muk.fits", uMap, overwrite=True)


