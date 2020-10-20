module PtsrcLens

using Distributed
using DelimitedFiles
using Random
using Setfield
using Statistics
using ProgressMeter
using HDF5
using FileIO
using DrWatson
using CUDA
using CMBLensing

export load_hdf5_cutouts, get_foreground_whitenoise_level, get_MAPs, sims, fluxcuts, get_fiducial_Cℓ, main_MAP_grid

const sims = 1:40
const fluxcuts = (2, 5, 10, Inf)


"""
Load the HDF5-converted cutouts (run make_h5output.jl to convert the
original text files to HDF5)
"""
function load_hdf5_cutouts()

    κs = @showprogress map(sims) do i
        ud_grade(FlatMap(h5read("data/sehgal_maps_h5/radio_sources/cutouts/kappa_sehgal_patch$i.h5","map"), θpix=1//2), 2)
    end
    
    ϕs = map(κs) do κ
        -2*(∇²\κ)
    end

    gs_radio = Dict(map([90,148]) do freq
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

    (;ϕs, κs, gs_ir, gs_radio, Ms_radio)

end


"""
Fit the mean white-noise level for radio and IR sources in units of μKarcmin.
"""
function get_foreground_whitenoise_level(;Ms_radio, gs_radio, gs_ir)

    μKarcmin_gs_radio = sort(Dict(map(collect(Ms_radio)) do ((survey,freq,fluxcut), Ms)
        (survey,freq,fluxcut) => if !(nothing in Ms)
            mean(sims) do i
                sqrt(mean(get_Cℓ(Ms[i] * gs_radio[freq][i], which=:QQ)[1000:2000]) / deg2rad(1/60)^2)
            end
        else
            missing
        end
    end))

    μKarcmin_gs_ir₀ = mean(sims) do i
        sqrt(mean(get_Cℓ(gs_ir[148][i], which=:QQ)[1000:2000]) / deg2rad(1/60)^2)
    end
    μKarcmin_gs_ir = Dict(
        (:deep,148,Inf) => μKarcmin_gs_ir₀,
        (:wide,148,Inf) => μKarcmin_gs_ir₀
    )

    (;μKarcmin_gs_radio, μKarcmin_gs_ir)

end


"""
Load fiducial Cℓs. For Cℓϕϕ, uses CAMB result up to ℓ=6000 and beyond
that uses the mean of the ϕ sims themselves (which have some
simulation-specific extra non-linear power).
"""
function get_fiducial_Cℓ(ϕs)

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
function get_MAPs(;Cℓ, Ms, gs, ϕs, noise_kwargs, ℓmax_data, polfrac_scale, ℓedges, μKarcmin_g, Nbatch=16, MAPs=nothing, sims=sims)
        
    @unpack ds = load_sim_dataset(;
        Cℓ = Cℓ,
        θpix = 2,
        Nside = 300,
        pol = :P,
        bandpass_mask = LowPass(ℓmax_data),
        Nbatch = Nbatch,
        noise_kwargs...
    )
    @unpack B = ds
    ds.Cϕ = Cℓ_to_Cov(ϕs[1], (Cℓ.total.ϕϕ, ℓedges, :Aϕ));

    Cℓg = noiseCℓs(μKarcminT=polfrac_scale*μKarcmin_g/√2, beamFWHM=0, ℓknee=0)
    Cg = Cℓ_to_Cov(Flat(Nside=300, θpix=2), Float32, S2, Cℓg.EE, Cℓg.BB)
    
    @showprogress pmap(sims) do sim

        sim′ = mod(sim,maximum(sims))+1

        Dict(map([

            (:nofg,     :fgcov, nothing,  nothing,  ϕs[sim]),
            (:corrfg,   :fgcov, gs[sim] , Ms[sim],  ϕs[sim]),
            (:uncorrfg, :fgcov, gs[sim′], Ms[sim′], ϕs[sim]),
            (:gaussfg,  :fgcov, nothing , 1,        ϕs[sim]),

            # (:nofg,     :nocov, nothing,  nothing,  ϕs[sim]),
            # (:corrfg,   :nocov, gs[sim] , Ms[sim],  ϕs[sim]),
            # (:uncorrfg, :nocov, gs[sim′], Ms[sim′], ϕs[sim]),
            # (:gaussfg,  :nocov, nothing,  1,        ϕs[sim]),
            
        ]) do (g_in_data, g_in_cov, g, M, ϕ)

            ds′ = resimulate(ds, ϕ=ϕ, seed=sim).ds
            if g_in_cov == :fgcov
                ds′.Cn += polfrac_scale^2*B*Cg*B'
            end
            if g_in_data in [:corrfg, :uncorrfg]
                ds′.d += polfrac_scale*B*M*g
            elseif g_in_data == :gaussfg
                g = simulate(Cg,seed=sim)
                ds′.d += polfrac_scale*B*M*g
            end
            ds′ = cu(ds′)

            (g_in_data,g_in_cov) => try
                if MAPs == nothing
                    fJ,ϕJ = MAP_joint(
                        ds′,
                        Nϕ       = :qe,
                        nsteps   = 30,
                        progress = false,
                    )
                else
                    fJ, ϕJ = cu.(MAPs[sim][g_in_data,g_in_cov][1:2])
                end
                gAϕ = gradient(Aϕ -> lnP(0,fJ,ϕJ,(Aϕ=Aϕ,),ds′), ones(Float32,length(ℓedges)-1))[1]
                (fJ=cpu(fJ), ϕJ=cpu(ϕJ), gAϕ)
            catch err
                rethrow(err)
                @warn "$sim $g_in_data $g_in_cov"
                nothing
            end

        end...)

    end
    
end



function main_MAP_grid(;surveys,freqs,ℓmax_datas,fluxcuts,polfrac_scales,Nbatch=16,overwrite=false,sims=sims)
    
    @unpack ϕs, κs, gs_ir, gs_radio, Ms_radio = load("data/sehgal_maps_h5/cutouts.jld2")
    @unpack (μKarcmin_gs_radio, μKarcmin_gs_ir) = get_foreground_whitenoise_level(;Ms_radio, gs_radio, gs_ir);
    Cℓ = get_fiducial_Cℓ(ϕs)

    ℓedges = [2:100:500; round.(Int, 10 .^ range(log10(502), log10(6000), length=10))]
    noises = Dict(
        (:deep, 90)  => (μKarcminT=0.68/√2, beamFWHM=2.3, ℓknee=200, αknee=2),
        (:wide, 90)  => (μKarcminT=2.9/√2,  beamFWHM=2.2, ℓknee=700, αknee=1.4),
        (:deep, 148) => (μKarcminT=0.96/√2, beamFWHM=1.5, ℓknee=200, αknee=2),
        (:wide, 148) => (μKarcminT=2.8/√2,  beamFWHM=1.4, ℓknee=700, αknee=1.4),
    )

    configs = collect(skipmissing(map(Iterators.product(surveys,freqs,ℓmax_datas,fluxcuts,polfrac_scales)) do (survey,freq,ℓmax_data,fluxcut,polfrac_scale)
        filename = datadir("MAPs", savename((;survey,freq,ℓmax_data,fluxcut,polfrac_scale,Nbatch), "jld2"))
        if overwrite || !isfile(filename)
            (survey,freq,ℓmax_data,fluxcut,polfrac_scale,filename)
        else
            missing
        end
    end))
    
    map(configs) do (survey,freq,ℓmax_data,fluxcut,polfrac_scale,filename)
        MAPs = get_MAPs(;
            Cℓ, Ms=Ms_radio[survey,freq,fluxcut], gs=gs_radio[freq], ϕs, noise_kwargs=noises[survey,freq], ℓmax_data, 
            polfrac_scale, ℓedges, μKarcmin_g=μKarcmin_gs_radio[survey,freq,fluxcut], Nbatch, sims
        )
        save(filename, "MAPs", MAPs)
    end

end


end