#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --gpus-per-task 1
#SBATCH -t 02:00:00
#SBATCH -C gpu
#SBATCH -A mp107
#SBATCH -o log/%x-%j.out
#=
srun -n 8 julia $(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}') "$ARGS"
exit
# =#


if length(ARGS) == 0

    using DrWatson

    sources = [:ir]
    surveys = [:deep, :wide]
    freqs = [148]
    ℓmax_datas = [3000, 4000, 5000]
    fluxcuts = [Inf]
    polfrac_scales = [1]
    Nbatch = 16
    overwrite = false

    configs = collect(skipmissing(map(Iterators.product(sources,surveys,freqs,ℓmax_datas,fluxcuts,polfrac_scales)) do (source,survey,freq,ℓmax_data,fluxcut,polfrac_scale)
        args = (;source,survey,freq,ℓmax_data,fluxcut,polfrac_scale,Nbatch)
        filename = datadir("MAPs", savename(args, "jld2"))
        cmd = "ARGS=\"$(args)\" sbatch $(@__FILE__)"
        if overwrite || !isfile(filename)
            println(cmd)
        else
            missing
        end
    end))

else

    using CUDA, CMBLensing, PtsrcLens, MPIClusterManagers
    CUDA.allowscalar(false)
    CMBLensing.init_MPI_workers()
    PtsrcLens.main_MAP_grid(;eval(Meta.parse(ARGS[1]))...)

end