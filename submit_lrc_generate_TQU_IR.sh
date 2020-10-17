#!/bin/bash
#SBATCH -J generate_TQU_IR
#SBATCH -N 1
#SBATCH -p lr3 # lr6 has 32 cores like cori
#SBATCH -q lr_normal
##SBATCH -q lr_debug
#SBATCH -t 23:59:59  # hh:mm:ss
##SBATCH -t 00:05:00  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eschaan@lbl.gov
#SBATCH -o /global/scratch/eschaan/SehgalRadioPolarizedSources/log/generate_TQU_IR.out
#SBATCH -e /global/scratch/eschaan/SehgalRadioPolarizedSources/log/generate_TQU_IR.err


cd /global/scratch/eschaan/SehgalRadioPolarizedSources/
source ~/python_profile.sh

# generate the full sky TQU maps
#source ./generate_TQU_IR.sh

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
#python generate_TQU_IR_cutouts.py "IRgal_S_11"
## BLAST very high flux population
#python generate_TQU_IR_cutouts.py "IRBlastPop"




## Syntax:
##python generate_TQU_IR_cutouts.py nu fluxCutmJy lKnee aKnee beamFwhm noiseT
#
## wide survey
## 148 GHz
#python generate_TQU_IR_cutouts.py 148.e9 10. 700. 1.4 1.4 2.
#python generate_TQU_IR_cutouts.py 148.e9 5. 700. 1.4 1.4 2.
#python generate_TQU_IR_cutouts.py 148.e9 2. 700. 1.4 1.4 2.
#
## deep survey
## 148 GHz
#python generate_TQU_IR_cutouts.py 148.e9 10. 200. 2. 1.5 0.68
#python generate_TQU_IR_cutouts.py 148.e9 5. 200. 2. 1.5 0.68
#python generate_TQU_IR_cutouts.py 148.e9 2. 200. 2. 1.5 0.68
