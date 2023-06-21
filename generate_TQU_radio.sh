#!/bin/bash


#####################################################################
# Fixed polarization fraction

# generate healpix map over the quadrant
#python generate_TQU_radio_hp.py 90.e9
#python generate_TQU_radio_hp.py 148.e9

# extract all the cutouts
#python generate_TQU_radio_cutouts.py 90.e9
#python generate_TQU_radio_cutouts.py 148.e9



# Generate all the point source masks
# Syntax
#python generate_TQU_radio_masks.py nu fluxCutmJy lKnee aKnee beamFwhm noiseT

# wide survey
# 90 GHz
#python generate_TQU_radio_masks.py 90.e9 10. 700. 1.4 2.2 2.
#python generate_TQU_radio_masks.py 90.e9 5. 700. 1.4 2.2 2.
python generate_TQU_radio_masks.py 90.e9 2. 700. 1.4 2.2 2.

# deep survey
# 90 GHz
#python generate_TQU_radio_masks.py 90.e9 10. 200. 2. 2.3 0.48
#python generate_TQU_radio_masks.py 90.e9 5. 200. 2. 2.3 0.48
#python generate_TQU_radio_masks.py 90.e9 2. 200. 2. 2.3 0.48


## wide survey
## 148 GHz
#python generate_TQU_radio_masks.py 148.e9 10. 700. 1.4 1.4 2.
#python generate_TQU_radio_masks.py 148.e9 5. 700. 1.4 1.4 2.
#python generate_TQU_radio_masks.py 148.e9 2. 700. 1.4 1.4 2.
#
## deep survey
## 148 GHz
#python generate_TQU_radio_masks.py 148.e9 10. 200. 2. 1.5 0.68
#python generate_TQU_radio_masks.py 148.e9 5. 200. 2. 1.5 0.68
#python generate_TQU_radio_masks.py 148.e9 2. 200. 2. 1.5 0.68

#####################################################################
# Gaussian polarization fraction

#python generate_TQU_radio_scatter_hp.py
