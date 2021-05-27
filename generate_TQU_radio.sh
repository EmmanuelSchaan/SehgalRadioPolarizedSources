#!/bin/bash


#####################################################################
# Fixed polarization fraction

python generate_TQU_radio_hp.py


# Syntax
#python generate_TQU_IR_cutouts.py nu fluxCutmJy lKnee aKnee beamFwhm noiseT

# wide survey
# 148 GHz
python generate_TQU_radio_cutouts.py 148.e9 10. 700. 1.4 1.4 2.
python generate_TQU_radio_cutouts.py 148.e9 5. 700. 1.4 1.4 2.
python generate_TQU_radio_cutouts.py 148.e9 2. 700. 1.4 1.4 2.

# deep survey
# 148 GHz
python generate_TQU_radio_cutouts.py 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_radio_cutouts.py 148.e9 5. 200. 2. 1.5 0.68
python generate_TQU_radio_cutouts.py 148.e9 2. 200. 2. 1.5 0.68


#####################################################################
# Gaussian polarization fraction

python generate_TQU_radio_scatter_hp.py


# Syntax
#python generate_TQU_IR_cutouts.py nu fluxCutmJy lKnee aKnee beamFwhm noiseT

# wide survey
# 148 GHz
python generate_TQU_radio_cutouts.py 148.e9 10. 700. 1.4 1.4 2.
python generate_TQU_radio_cutouts.py 148.e9 5. 700. 1.4 1.4 2.
python generate_TQU_radio_cutouts.py 148.e9 2. 700. 1.4 1.4 2.

# deep survey
# 148 GHz
python generate_TQU_radio_cutouts.py 148.e9 10. 200. 2. 1.5 0.68
python generate_TQU_radio_cutouts.py 148.e9 5. 200. 2. 1.5 0.68
python generate_TQU_radio_cutouts.py 148.e9 2. 200. 2. 1.5 0.68

