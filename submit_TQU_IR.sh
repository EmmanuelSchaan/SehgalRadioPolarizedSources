#!/bin/bash
#SBATCH -J TQU_IR
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use haswell/knl nodes
#SBATCH -t 11:59:59  #15:59:59  # hh:mm:ss
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=eschaan@lbl.gov
#SBATCH -o /global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/log/tqu_ir_hist.out
#SBATCH -e /global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/log/tqu_ir_hist.err


cd /global/cscratch1/sd/eschaan/SehgalRadioPolarizedSources/
source ~/python_profile.sh

source ./generate_TQU_IR.sh

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


