
using DrWatson

sources = [:radio]
surveys = [:deep, :wide]
freqs = [90]
ℓmax_datas = [5000]
fluxcuts = [5]
polfrac_scales = [1]
Nbatch = 16
overwrite = false

configs = collect(skipmissing(map(Iterators.product(sources,surveys,freqs,ℓmax_datas,fluxcuts,polfrac_scales)) do (source,survey,freq,ℓmax_data,fluxcut,polfrac_scale)
    args = (;source,survey,freq,ℓmax_data,fluxcut,polfrac_scale,Nbatch)
    filename = datadir("MAPs", savename(args, "jld2"))
    cmd = "sbatch -C gpu -t 04:00:00 -c 2 -G 4 -n 5 -N 1 -A mp107 -o log/%j.out --wrap=\"srun julia -e 'using CUDA, CMBLensing, PtsrcLens, MPIClusterManagers; CMBLensing.init_MPI_workers(); PtsrcLens.main_MAP_grid(;$(args)...)'\""
    if overwrite || !isfile(filename)
        println(cmd)
    else
        missing
    end
end))
