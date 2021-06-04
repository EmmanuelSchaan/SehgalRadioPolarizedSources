#!/usr/bin/env sh
#SBATCH -G 4
#SBATCH -n 5
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -C gpu
#SBATCH -A m1759
#SBATCH -o log/%x-%j.out
#=
srun julia $(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
exit
# =#


using CUDA, CMBLensing, PtsrcLens, MPIClusterManagers
CMBLensing.init_MPI_workers()
PtsrcLens.main_MAP_grid(;
    # surveys = [:deep, :wide],
    surveys = [:wide],
    freqs = [90, 148],
    ℓmax_datas = [3000, 5000],
    fluxcuts = [2, 5, 10],
    polfrac_scales = [1],
    Nbatch = 16,
    overwrite = false
)