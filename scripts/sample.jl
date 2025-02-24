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
using NCDatasets

s = ArgParseSettings()

@add_arg_table s begin
    "--nmcmc"
        help = "Number of MCMC iterations"
        default = 1000
        arg_type = Int
    "--npost"
        help = "Number of posterior samples to draw"
        default = 256
        arg_type = Int
    "--nsel"
        help = "Number of detected injections to use to estimate the selection normalization"
        default = 4096
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

    "--o4a"
	help = "Whether or not to include GW230529"
        default = false
        arg_type = Bool

    "--mlow"
        help = "Minimum mass for event selection"
	default = nothing
	arg_type = Float64
end

parsed_args = parse_args(s)

Nmcmc = parsed_args["nmcmc"]
Nchain = Threads.nthreads()
Npost = parsed_args["npost"]
Nsel = parsed_args["nsel"]
target_accept = parsed_args["target_accept"]
m_low = parsed_args["mlow"]
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

log_dN_default = make_log_dN(PowerLawGaussian(), PowerLawPairing(), 3.2, -2.2, 9.8, 1.6, 0.77, -1.7, 0.14, 1.1, 0.94; m_low = m_low)
a = log_dN_default(10,3)

sampso3a, fnameso3a, gwnameso3a = load_pe_samples("../../o3a", "*GW*_nocosmo.h5", "C01:Mixed/posterior_samples")
sampso3b, fnameso3b, gwnameso3b = load_pe_samples("../../o3b", "*GW*_nocosmo.h5", "C01:Mixed/posterior_samples")
samps190425, fnames190425, gwnames190425 = load_pe_samples("../../o3a", "*GW*_nocosmo.h5", "C01:IMRPhenomPv2_NRTidal:HighSpin/posterior_samples")
println(gwnames190425)

if parsed_args["o4a"]
  gwtc_table = DataFrame(CSV.File("../../GWTC-up.csv"))

  sampso4a, fnameso4a, gwnameso4a = load_pe_samples("../../o4a", "*GW*_nocosmo.h5", "Combined_PHM_highSpin/posterior_samples")
  
  samps = vcat(sampso3a, sampso3b, samps190425, sampso4a)
  fnames = vcat(fnameso3a, fnameso3b, fnames190425, fnameso4a)
  gwnames = vcat(gwnameso3a, gwnameso3b, gwnames190425, gwnameso4a)
else
  gwtc_table = DataFrame(CSV.File("../../GWTC.csv"))


  samps = vcat(sampso3a, sampso3b, samps190425)
  fnames = vcat(fnameso3a, fnameso3b, fnames190425)
  gwnames = vcat(gwnameso3a, gwnameso3b, gwnames190425)
end

narrow_samps, narrow_fnames, narrow_gwnames = filter_selected(samps, fnames, gwnames, gwtc_table; m_low = m_low)

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

m1s_sel, m2s_sel, zs_sel, pdraw, Ndraw = read_selection("../../o1+o2+o3_mixture_real+semianalytic-LIGO-T2100377-v2.hdf5")
pdraw = pdraw ./ (md_sfr_zwt.(zs_sel) ./ md_sfr(zref))
m1s_sel, m2s_sel, pdraw, Ndraw = resample_selection(log_dN_default, m1s_sel, m2s_sel, pdraw, Ndraw)
@info "Number of resampled selection function samples = $(Ndraw)"

f = Nsel / length(m1s_sel)
Ndraw = round(Int, f*Ndraw)
m1s_sel = m1s_sel[1:Nsel]
m2s_sel = m2s_sel[1:Nsel]
pdraw = pdraw[1:Nsel]

model = model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log.(pdraw), Ndraw, m_low = m_low)

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

if parsed_args["o4a"]
  suffix = suffix*"_including_230529"
end
if m_low !== nothing
  suffix = string(suffix*"_", round(m_low; digits = 3))
end
to_netcdf(idata, joinpath(@__DIR__, "..", "chains", "chain" * suffix * ".nc"))
