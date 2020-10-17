using Distributed
addprocs(8)

##

@everywhere begin
    using HDF5
    using Glob
    using DelimitedFiles
    using ProgressMeter
end

##

dirs = [
    "dat/sehgal_maps/radio_sources/cutouts/*.txt",
    "dat/sehgal_maps/ir_sources/cutouts/*.txt",
]

@showprogress pmap(mapreduce(glob, vcat, dirs)) do f_in
    f_out = replace(splitext(f_in)[1]*".h5", "sehgal_maps"=>"sehgal_maps_h5")
    if !isfile(f_out)
        mkpath(dirname(f_out))
        h5write(f_out, "map", Float32.(readdlm(f_in, use_mmap=false)))
    end
end;


