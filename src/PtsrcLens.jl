module PtsrcLens

using Adapt
using CMBLensing
using CUDA
using DataStructures
using DelimitedFiles
using Distributed
using DrWatson
using FileIO
using HDF5
using Match
using Measurements: value, uncertainty, ±
using Memoization
using ProgressMeter
using Random
using Setfield
using Statistics

export load_hdf5_cutouts, load_jld2_cutouts, get_foreground_noise, get_MAPs, sims, fluxcuts, 
    get_fiducial_Cℓ, main_MAP_grid, noises, FgDataContrib, load_ptsrclens_dataset

    
const sims = 1:40
const fluxcuts = (2, 5, 10, Inf)
const noises = Dict(
    (:deep, 90)  => (μKarcminT=0.68/√2, beamFWHM=2.3, ℓknee=200, αknee=2),
    (:wide, 90)  => (μKarcminT=2.9/√2,  beamFWHM=2.2, ℓknee=700, αknee=1.4),
    (:deep, 148) => (μKarcminT=0.96/√2, beamFWHM=1.5, ℓknee=200, αknee=2),
    (:wide, 148) => (μKarcminT=2.8/√2,  beamFWHM=1.4, ℓknee=700, αknee=1.4),
)
const fskys = Dict(
    :deep => 0.03, 
    :wide => 0.6
)
const DefaultRNG = Xoshiro
  
"""
Load the HDF5-converted cutouts (run make_h5output.jl to convert the
original text files to HDF5)
"""
function load_hdf5_cutouts(;freqs_radio=[90,148])

    κs = @showprogress map(sims) do i
        ud_grade(FlatMap(h5read("data/sehgal_maps_h5/radio_sources/cutouts/kappa_sehgal_patch$i.h5","map"), θpix=1//2), 2)
    end
    
    ϕs = map(κs) do κ
        -2*(∇²\κ)
    end

    gs_radio = Dict(map(freqs_radio) do freq
        freq => @showprogress map(sims) do i 
            ud_grade(
                FlatQUMap((h5read("data/sehgal_maps_h5/radio_sources/cutouts/ps_sehgal_$(freq)ghz_$(pol)_patch$i.h5", "map") for pol in "QU")..., θpix=1//2), 
                2, 
                anti_aliasing=false, 
                deconv_pixwin=false
            )
        end
    end)

    gs_ir = Dict(148 => @showprogress map(sims) do sim
        mapreduce(+, ["gal_S_$(S)" for S=1:9]) do bin
            ud_grade(
                FlatQUMap((Float32.(h5read("data/sehgal_maps_h5/ir_sources/cutouts/ir_IR$(bin)_sehgal_$(pol)_patch$(sim).h5", "map")) for pol="QU")..., θpix=1//2), 2,
                anti_aliasing=false, deconv_pixwin=false
            )
        end
    end)

    surveykeys = Dict(
        (90,  :deep) => "beam2.3_noise0.48_lknee200_aknee2.0",
        (90,  :wide) => "beam2.2_noise2.0_lknee700_aknee1.4",
        (148, :deep) => "beam1.5_noise0.68_lknee200_aknee2.0",
        (148, :wide) => "beam1.4_noise2.0_lknee700_aknee1.4"
    )

    Ms_radio = Dict(map(Iterators.product([:deep,:wide],[90,148],fluxcuts)) do (survey,freq,fluxcut)
        (survey,freq,fluxcut) => if fluxcut == Inf
            one.(sims)
        else
            @showprogress map(sims) do i
                filename = "data/sehgal_maps_h5/radio_sources/cutouts/ps_mask_$(freq)ghz_$(fluxcut)mJy_$(surveykeys[freq,survey])_T_patch$(i).h5"
                if isfile(filename)
                    M = ud_grade(FlatMap(Float32.(h5read(filename,"map")),θpix=1//2), 2, anti_aliasing=false, deconv_pixwin=false)
                    Diagonal(FlatQUMap(M,M))
                else
                    nothing
                end
            end
        end
    end)

    gs = Dict(:radio => gs_radio, :ir => gs_ir)
    Ms = Dict(:radio => Ms_radio, :ir => DefaultDict(DefaultDict(1)));

    (;ϕs, κs, gs, Ms)

end

"""
Load the JLD2-converted cutouts (output from load_hdf5_cutouts)
"""
@memoize load_jld2_cutouts(filename="data/sehgal_maps_h5/cutouts.jld2") = load(filename) 


"""
Fit the mean white-noise level for radio and IR sources in units of μKarcmin.
"""
@memoize function get_foreground_noise(;Ms, gs, ℓrange=5000:10000)

    fg_noise_radio = sort(Dict(map(collect(Ms[:radio])) do ((survey,freq,fluxcut), Ms)
        
        (survey, freq, fluxcut) => if haskey(gs[:radio], freq)

            fg_noise = mean(sims) do i
                sqrt(mean(get_Cℓ(Ms[i] * gs[:radio][freq][i], which=:QQ)[1000:2000]) / deg2rad(1/60)^2)
            end

            fsky = fskys[survey]
            Bℓ = beamCℓs(;noises[survey,freq].beamFWHM,ℓmax=last(ℓrange))[ℓrange]
            Cℓfg = fg_noise^2
            Cℓnoise = noises[survey,freq].μKarcminT^2

            F = (fsky/2) * sum(@. (2ℓrange+1) * (Bℓ * Cℓfg)^2 / (Bℓ * Cℓfg + Cℓnoise)^2)
            
            fg_noise ± fg_noise * sqrt(inv(F))/2

        end

    end))

    fg_noise_ir₀ = mean(sims) do i
        sqrt(mean(get_Cℓ(gs[:ir][148][i], which=:QQ)[1000:2000]) / deg2rad(1/60)^2)
    end
    fg_noise_ir = Dict(
        (:deep,148,Inf) => fg_noise_ir₀,
        (:wide,148,Inf) => fg_noise_ir₀
    )

    Dict(:radio => fg_noise_radio, :ir => fg_noise_ir)

end


"""
Load fiducial Cℓs. For Cℓϕϕ, uses CAMB result up to ℓ=6000 and beyond
that uses the mean of the ϕ sims themselves (which have some
simulation-specific extra non-linear power).
"""
@memoize function get_fiducial_Cℓ(ϕs)

    Cℓ = camb(r=0.001, ωb=0.02268, ωc=0.1081, nₛ=0.961, H0=72.4, θs=nothing, logA=log(2.41*10), k_pivot=0.05, ℓmax=10000)

    LP = LowPass(1100,Δℓ=200).diag.Wℓ
    HP = HighPass(901,Δℓ=200).diag.Wℓ
    Cℓϕϕ = CMBLensing.extrapolate_Cℓs(
        2:10000, 
        2:6000, 
        nan2zero.((LP * Cℓ.total.ϕϕ)[2:6000]) + nan2zero.((HP * CMBLensing.smooth(mean(get_Cℓ(ϕ) for ϕ in ϕs), xscale=:log, yscale=:log))[2:6000])
    )

    @set! Cℓ.total.ϕϕ = Cℓϕϕ

    Cℓ

end


"""
Compute the joint MAP for some sims.
"""
function get_MAPs(;
    Cℓ, Ms, gs, ϕs, noise_kwargs, ℓmax_data, polfrac_scale, ℓedges, fg_noise,
    Nbatch=16, MAPs=nothing, sims=sims
)

    pbar = Progress(5*length(sims), 1, "get_MAPs: ")
    ProgressMeter.update!(pbar)
    update_pbar = RemoteChannel(()->Channel{Bool}(), 1)
    @async while take!(update_pbar)
        next!(pbar)
    end

    @unpack ds, proj = load_sim(;
        Cℓ = Cℓ,
        θpix = 2,
        Nside = 300,
        pol = :P,
        bandpass_mask = LowPass(ℓmax_data),
        Nbatch = Nbatch,
        noise_kwargs...
    )
    @unpack B = ds
    T = real(eltype(ds.d))
    ds.Cϕ = Cℓ_to_Cov(:I, proj, (Cℓ.total.ϕϕ, ℓedges, :Aϕ))

    Cℓg = noiseCℓs(μKarcminT=polfrac_scale*value(fg_noise)/√2, beamFWHM=0, ℓknee=0)
    Cg = Cℓ_to_Cov(:P, proj, Cℓg.EE, Cℓg.BB)

    pool = CachingPool(procs())
    
    pmap(pool, sims) do sim

        sim′ = mod(sim,maximum(sims))+1

        Dict(map([

            (:nofg,      :nofgcov,  nothing,  nothing,  ϕs[sim], T(polfrac_scale)),
            (:nofg,      :fgcov,    nothing,  nothing,  ϕs[sim], T(polfrac_scale)),
            (:corrfg,    :fgcov,    gs[sim] , Ms[sim],  ϕs[sim], T(polfrac_scale)),
            (:uncorrfg,  :fgcov,    gs[sim′], Ms[sim′], ϕs[sim], T(polfrac_scale)),
            (:corrfg,    :nofgcov,  gs[sim] , Ms[sim],  ϕs[sim], T(polfrac_scale)),
            (:gaussfg,   :fgcov,    nothing , 1,        ϕs[sim], T(polfrac_scale)),
            (:gaussfg2,  :fgcov,    nothing , 1,        ϕs[sim], T(polfrac_scale * (1 + uncertainty(fg_noise)/value(fg_noise)))),

        ]) do (g_in_data, g_in_cov, g, M, ϕ, polfrac_scale)

            ds′ = copy(ds)
            ds′.d = simulate(DefaultRNG(sim), ds, ϕ=ϕ).d
            if g_in_cov == :fgcov
                ds′.Cn += polfrac_scale^2*B*Cg*B'
            end
            if g_in_data in [:corrfg, :uncorrfg]
                ds′.d += polfrac_scale*B*M*g
            elseif g_in_data in [:gaussfg, :gaussfg2]
                g = simulate(DefaultRNG(sim), Cg)
                ds′.d += polfrac_scale*B*M*g
            end
            ds′ = cu(ds′)

            (g_in_data, g_in_cov) => try
                if MAPs == nothing
                    fJ,ϕJ = MAP_joint(
                        ds′,
                        nsteps   = 30,
                        progress = false,
                        αtol     = 1e-6, 
                        nburnin_update_hessian = Inf
                    )
                else
                    fJ, ϕJ = cu.(MAPs[sim][g_in_data,g_in_cov][1:2])
                end
                gAϕ = gradient(Aϕ -> lnP(0,fJ,ϕJ,(Aϕ=Aϕ,),ds′), ones(Float32,length(ℓedges)-1))[1]
                put!(update_pbar, true)
                (fJ=cpu(fJ), ϕJ=cpu(ϕJ), gAϕ)
            catch err
                if (err isa InterruptException) || (err isa RemoteException && err.captured isa CapturedException && err.captured.ex isa InterruptException)
                    rethrow(err)
                else
                    @warn "$sim $g_in_data $g_in_cov" err
                    nothing
                end
            end

        end...)

    end
    
end



function main_MAP_grid(;source, survey, freq, ℓmax_data, fluxcut, polfrac_scale, Nbatch=16, overwrite=false, sims=sims)
    
    @show (source, survey, freq, ℓmax_data, fluxcut, polfrac_scale)

    @unpack ϕs, κs, gs, Ms = load("data/sehgal_maps_h5/cutouts.jld2")
    fg_noise = get_foreground_noise(;Ms, gs);
    Cℓ = get_fiducial_Cℓ(ϕs)

    ℓedges = [2; 100:50:500; round.(Int, 10 .^ range(log10(502), log10(5000), length=40))];

    MAPs = get_MAPs(;
        Cℓ, Ms=Ms[source][survey,freq,fluxcut], gs=gs[source][freq], ϕs, noise_kwargs=noises[survey,freq], ℓmax_data, 
        polfrac_scale, ℓedges, fg_noise=fg_noise[source][survey,freq,fluxcut], Nbatch, sims
    )

    filename = datadir("MAPs", savename((;source,survey,freq,ℓmax_data,fluxcut,polfrac_scale,Nbatch), "jld2"))
    @time save(filename, "MAPs", MAPs)

end



@enum FgDataContrib CUTOUT_CORRELATED CUTOUT_UNCORRELATED GAUSSIAN GAUSSIAN_OFF_BY_SIGMA

function load_ptsrclens_dataset(;
    source,
    survey,
    freq,
    ℓmax_data,
    fluxcut,
    polfrac_scale,
    sim,
    fg_data_contrib :: Union{Nothing,FgDataContrib},
    fg_cov_modeled :: Bool,
    Nbatch = 16,
    ℓedges_ϕ = [2; 100:50:500; round.(Int, 10 .^ range(log10(502), log10(5000), length=40))],
    ℓedges_E = 500:100:5000,
    T = Float32,
    storage = CUDA.functional() ? CuArray : Array,
)

    @unpack ϕs, κs, gs, Ms = load_jld2_cutouts()
    Cℓ = get_fiducial_Cℓ(ϕs)
    fg_noises = get_foreground_noise(;Ms, gs)

    noise_kwargs = noises[survey,freq]
    fg_noise = fg_noises[source][survey,freq,fluxcut]

    @unpack ds, proj = load_sim(;
        Cℓ = Cℓ,
        θpix = 2,
        Nside = 300,
        pol = :P,
        bandpass_mask = LowPass(ℓmax_data),
        Nbatch,
        T,
        noise_kwargs...
    )
    @unpack B,M = ds
    T = real(eltype(ds.d))
    ds.Cϕ = Cℓ_to_Cov(:I, proj, (Cℓ.total.ϕϕ, ℓedges_ϕ, :Aϕ))

    Cℓg = noiseCℓs(μKarcminT=polfrac_scale*value(fg_noise)/√2, beamFWHM=0, ℓknee=0)
    Cg = Cℓ_to_Cov(:P, proj, Cℓg.EE, Cℓg.BB)

    ϕ = ϕs[sim]
    sim′ = mod(sim,maximum(sims))+1
    g = @match string(fg_data_contrib) begin
        "CUTOUT_CORRELATED"     => gs[source][freq][sim]
        "CUTOUT_UNCORRELATED"   => gs[source][freq][sim′]
        "GAUSSIAN"              => simulate(DefaultRNG(sim),Cg)
        "GAUSSIAN_OFF_BY_SIGMA" => simulate(DefaultRNG(sim),Cg) * (1 + uncertainty(fg_noise)/value(fg_noise))
    end
    Mg = @match string(fg_data_contrib) begin
        "CUTOUT_CORRELATED"   => Ms[source][survey,freq,fluxcut][sim]
        "CUTOUT_UNCORRELATED" => Ms[source][survey,freq,fluxcut][sim′]
        _                     => 1
    end

    ds.d = simulate(DefaultRNG(sim), ds, ϕ=ϕ).d
    if g != nothing 
        ds.d += polfrac_scale*M*B*Mg*g
    end
    if fg_cov_modeled
        ds.Cn += polfrac_scale^2*M*B*Cg*B'*M'
    end

    adapt(storage, (; ds, ϕ, g, Mg))

end



end