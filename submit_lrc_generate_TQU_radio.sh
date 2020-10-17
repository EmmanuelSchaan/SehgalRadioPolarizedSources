#!/bin/bash
#SBATCH -J generate_TQU_radio
#SBATCH -N 1
#SBATCH -p lr3 # lr6 has 32 cores like cori
#SBATCH -q lr_normal
##SBATCH -q lr_debug
#SBATCH -t 23:59:59  # hh:mm:ss
##SBATCH -t 00:05:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eschaan@lbl.gov
#SBATCH -o /global/scratch/eschaan/SehgalRadioPolarizedSources/log/generate_TQU_radio.out
#SBATCH -e /global/scratch/eschaan/SehgalRadioPolarizedSources/log/generate_TQU_radio.err


cd /global/scratch/eschaan/SehgalRadioPolarizedSources/

# Syntax:
#python generate_TQU_radio_cutouts.py nu fluxCutmJy lKnee aKnee beamFwhm noiseT

# wide survey
# 90 GHz
#python generate_TQU_radio_cutouts.py 90.e9 10. 700. 1.4 2.2 2.
#python generate_TQU_radio_cutouts.py 90.e9 5. 700. 1.4 2.2 2.
#python generate_TQU_radio_cutouts.py 90.e9 2. 700. 1.4 2.2 2.
# 148 GHz
#python generate_TQU_radio_cutouts.py 148.e9 10. 700. 1.4 1.4 2.
#python generate_TQU_radio_cutouts.py 148.e9 5. 700. 1.4 1.4 2.
#python generate_TQU_radio_cutouts.py 148.e9 2. 700. 1.4 1.4 2.

# deep survey
# 90 GHz
python generate_TQU_radio_cutouts.py 90.e9 10. 200. 2. 2.3 0.48
python generate_TQU_radio_cutouts.py 90.e9 5. 200. 2. 2.3 0.48
python generate_TQU_radio_cutouts.py 90.e9 2. 200. 2. 2.3 0.48
# 148 GHz
#python generate_TQU_radio_cutouts.py 148.e9 10. 200. 2. 1.5 0.68
#python generate_TQU_radio_cutouts.py 148.e9 5. 200. 2. 1.5 0.68
#python generate_TQU_radio_cutouts.py 148.e9 2. 200. 2. 1.5 0.68
