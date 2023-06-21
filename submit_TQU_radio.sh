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
source ~/python3_profile.sh

source ./generate_TQU_radio.sh



