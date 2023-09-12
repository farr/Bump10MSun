using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ArgParse
using ArviZ
using Bump10MSun
using CSV
using DataFrames
using Glob
using HDF5
using Logging
using MCMCChainsStorage # Needed for `append_generated_quantities`
using Printf
using Random
using StatsBase
using Turing

s = ArgParseSettings()

@add_arg_table s begin
    "--nmcmc"
        help = "Number of MCMC iterations"
        default = 1000
        arg_type = Int
    "--npost"
        help = "Number of posterior samples to draw"
        default = 512
        arg_type = Int
    "--nsel"
        help = "Number of detected injections to use to estimate the selection normalization"
        default = 6500
        arg_type = Int
    "--model"
        help = "Model to fit, one of [broken_pl, two_broken_pl, broken_pl_plp, two_broken_pl_plp]"
        default = "broken_pl"
        arg_type = String
        required = true
    "--target_accept"
        help = "Target acceptance rate for NUTS"
        default = 0.85
        arg_type = Float64
    "--seed"
        help = "RNG seed"
        default = 8572587030709747370
        arg_type = Int
    "--extra-events"
        help = "Include by hand extra events"
        action = :store_true
end

parsed_args = parse_args(s)

Nmcmc = parsed_args["nmcmc"]
Nchain = Threads.nthreads()
Npost = parsed_args["npost"]
Nsel = parsed_args["nsel"]
target_accept = parsed_args["target_accept"]
Random.seed!(parsed_args["seed"])

model_functions = Dict(
    "broken_pl" => broken_pl_plp_model,
    "broken_pl_gaussian" => power_law_plus_gaussian_model
)
model_suffixes = Dict(
    "broken_pl" => "_bpl",
    "broken_pl_gaussian" => "_bplg"
)

if parsed_args["model"] in keys(model_functions)
    model = model_functions[parsed_args["model"]]
    suffix = model_suffixes[parsed_args["model"]]
else
    error("--model argument must be one of $(keys(model_functions)); got $(parsed_args["model"])")
end

if parsed_args["extra-events"]
    suffix = suffix * "_extra"
end

log_dN_default = make_log_dN(PowerLawGaussian(), PowerLawPairing(), -1.5, -0.85, 9.7, 1.1, 0.55, 4.4, 4.6, 1.1, 1.9)

sampso3a, fnameso3a, gwnameso3a = load_pe_samples("/Users/wfarr/Research/gwtc-2.1", "*GW*_nocosmo.h5", "C01:Mixed/posterior_samples")
sampso3b, fnameso3b, gwnameso3b = load_pe_samples("/Users/wfarr/Research/o3b_data/PE", "*GW*_nocosmo.h5", "C01:Mixed/posterior_samples")
gwtc3_table = DataFrame(CSV.File("/Users/wfarr/Research/gwtc-3-tables/GWTC.csv"))

samps = vcat(sampso3a, sampso3b)
fnames = vcat(fnameso3a, fnameso3b)
gwnames = vcat(gwnameso3a, gwnameso3b)

narrow_samps, narrow_fnames, narrow_gwnames = filter_selected(samps, fnames, gwnames, gwtc3_table)

if parsed_args["extra-events"]
    df = DataFrame(CSV.File("/Users/wfarr/Research/o4a_data/pe/S230529ay/EXP6_pesummary.dat.gz"))

    nt_array = [(mass_1_source=m1, mass_2_source=m2, redshift=z) for (m1, m2, z) in zip(df[:, :mass_1_source], df[:, :mass_2_source], df[:, :redshift])]

    push!(narrow_samps, nt_array)
    push!(narrow_fnames, "Exp6_pesummary.dat.gz")
    push!(narrow_gwnames, "S230529ay")
end

@info "Analyzing $(length(narrow_gwnames)) events:"
for n in narrow_gwnames
    @info "  $(n)"
end

pop_wts_and_wts = map(narrow_samps) do samps
    map(samps) do s
        m1 = s.mass_1_source
        m2 = s.mass_2_source
        z = s.redshift
        pop_wt = exp.(log_dN_default(m1, m2)).*md_sfr_zwt(z)
        (pop_wt, pop_wt ./ li_prior_wt(m1, m2, z))
    end
end
pop_wts = [[x[1] for x in xx] for xx in pop_wts_and_wts]
wts = [[x[2] for x in xx] for xx in pop_wts_and_wts]

neffs = [sum(w)^2 / sum(w.^2) for w in wts]
@info @sprintf("Minimum PE Neff = %.1f", minimum(neffs))
inds = [sample(1:length(w), weights(w), Npost) for w in wts]

@info @sprintf("Maximum of 99%% redshift in sample is %.1f", maximum([quantile([x.redshift for x in s], 0.99) for s in narrow_samps]))

m1s = [[xs[i].mass_1_source for i in is] for (xs, is) in zip(narrow_samps, inds)]
m2s = [[xs[i].mass_2_source for i in is] for (xs, is) in zip(narrow_samps, inds)]
log_wts = [[log(w[i]) for i in is] for (w, is) in zip(pop_wts, inds)]

m1s_sel, m2s_sel, zs_sel, pdraw, Ndraw = read_selection("/Users/wfarr/Research/o3b_data/O1O2O3-Sensitivity/o1+o2+o3_mixture_real+semianalytic-LIGO-T2100377-v2.hdf5")
pdraw = pdraw ./ (md_sfr_zwt.(zs_sel) ./ md_sfr(zref))
m1s_sel, m2s_sel, pdraw, Ndraw = resample_selection(log_dN_default, m1s_sel, m2s_sel, pdraw, Ndraw)
@info "Number of resampled selection function samples = $(Ndraw)"

f = Nsel / length(m1s_sel)
Ndraw = round(Int, f*Ndraw)
m1s_sel = m1s_sel[1:Nsel]
m2s_sel = m2s_sel[1:Nsel]
pdraw = pdraw[1:Nsel]

model = model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log.(pdraw), Ndraw)

if Nchain == 1
    trace = sample(model, NUTS(Nmcmc, target_accept), Nmcmc)
else
    trace = sample(model, NUTS(Nmcmc, target_accept), MCMCThreads(), Nmcmc, Nchain)
end
trace = with_logger(NullLogger()) do 
    append_generated_quantities(trace, generated_quantities(model, trace))
end
@info @sprintf("Minimum Neff_sel = %.1f, 4*nobs = %d", minimum(trace[:Neff_sel]), 4*length(narrow_gwnames))
@info @sprintf("Minimum Neff_samps = %.1f", minimum([minimum(trace[n]) for n in namesingroup(trace, :Neff_samps)]))

idata = from_mcmcchains(trace, coords=Dict(:gwnames => narrow_gwnames), dims=(Neff_samps=(:gwnames,), m1s_popwt=(:gwnames,), m2s_popwt=(:gwnames,)), library="Turing")

to_netcdf(idata, joinpath(@__DIR__, "..", "chains", "chain" * suffix * ".nc"))