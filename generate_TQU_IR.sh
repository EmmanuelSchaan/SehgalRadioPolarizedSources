#!/bin/bash
# WARNING:
# the python script called here will abort right away
# if none of the sources in the catalog are below 100mJy at 148 GHz.

# Basic population of IR galaxies
# the files are ranked in order of growing galaxy IR flux
# No need to do the high flux ones!
#
# I believe the files split the sources based on their log10 flux at 350GHz,
# so the flux bounds in each file should be:
#In [10]: np.logspace(np.log10(1.e-2), np.log10(3.e2), 12, 10.)
#Out[10]:
#array([1.00000000e-02, 2.55276610e-02, 6.51661474e-02, 1.66353932e-01,
#       4.24662676e-01, 1.08406448e+00, 2.76736306e+00, 7.06443058e+00,
#       1.80338389e+01, 4.60361725e+01, 1.17519580e+02, 3.00000000e+02])
# Actually, this doesn't seem to really be the case...


# First, generate the healpix maps
#python generate_TQU_IR_hp.py "IRgal_S_1"
#python generate_TQU_IR_hp.py "IRgal_S_2"
#python generate_TQU_IR_hp.py "IRgal_S_3"
#python generate_TQU_IR_hp.py "IRgal_S_4"
#python generate_TQU_IR_hp.py "IRgal_S_5"
#python generate_TQU_IR_hp.py "IRgal_S_6"
#python generate_TQU_IR_hp.py "IRgal_S_7"
#python generate_TQU_IR_hp.py "IRgal_S_8"
#python generate_TQU_IR_hp.py "IRgal_S_9"
#python generate_TQU_IR_hp.py "IRgal_S_10"
#python generate_TQU_IR_hp.py "IRgal_S_11" # got an error with this one
## BLAST very high flux population
#python generate_TQU_IR_hp.py "IRBlastPop"


# Then, extract the cutouts
#python generate_TQU_IR_cutouts.py "IRgal_S_1"
#python generate_TQU_IR_cutouts.py "IRgal_S_2"
#python generate_TQU_IR_cutouts.py "IRgal_S_3"
#python generate_TQU_IR_cutouts.py "IRgal_S_4"
#python generate_TQU_IR_cutouts.py "IRgal_S_5"
#python generate_TQU_IR_cutouts.py "IRgal_S_6"
#python generate_TQU_IR_cutouts.py "IRgal_S_7"
#python generate_TQU_IR_cutouts.py "IRgal_S_8"
#python generate_TQU_IR_cutouts.py "IRgal_S_9"
#python generate_TQU_IR_cutouts.py "IRgal_S_10"
#python generate_TQU_IR_cutouts.py "IRgal_S_11" # did not run this one
# BLAST very high flux population
#python generate_TQU_IR_cutouts.py "IRBlastPop"


# Then create the PS masks
# Syntax
#python generate_TQU_IR_masks.py sourceCatalog nu fluxCutmJy lKnee aKnee beamFwhm noiseT

##python generate_TQU_IR_masks.py "IRgal_S_1" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_2" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_3" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_4" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_5" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_6" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_7" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_8" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_9" 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_IR_masks.py "IRgal_S_10" 148.e9  10. 200. 2. 1.5 0.68
#python generate_TQU_IR_masks.py "IRgal_S_11" 148.e9 10. 200. 2. 1.5 0.68 # did not run this one
# BLAST very high flux population
python generate_TQU_IR_masks.py "IRBlastPop" 148.e9 10. 200. 2. 1.5 0.68

