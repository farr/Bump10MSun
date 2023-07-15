using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ArgParse
using ArviZ
using Bump10MSun
using Glob
using HDF5
using Logging
using MCMCChainsStorage
using Printf
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
        default = 256
        arg_type = Int
    "--nsel"
        help = "Number of samples to select"
        default = 4096
        arg_type = Int
    "--model"
        help = "Model to fit, one of [broken_pl, two_broken_pl, broken_pl_plp, two_broken_pl_plp]"
        default = "broken_pl"
        arg_type = String
    "--target_accept"
        help = "Target acceptance rate for HMC"
        default = 0.85
        arg_type = Float64
end

parsed_args = parse_args(s)

Nmcmc = parsed_args["nmcmc"]
Nchain = Threads.nthreads()
Npost = parsed_args["npost"]
Nsel = parsed_args["nsel"]
target_accept = parsed_args["target_accept"]

model_functions = Dict(
    "broken_pl" => broken_pl_model,
    "two_broken_pl" => two_broken_pl_model,
    "broken_pl_plp" => broken_pl_plp_model,
    "two_broken_pl_plp" => two_broken_pl_plp_model
)
model_suffixes = Dict(
    "broken_pl" => "",
    "two_broken_pl" => "_tb",
    "broken_pl_plp" => "_plp",
    "two_broken_pl_plp" => "_tb_plp"
)

if parsed_args["model"] in keys(model_functions)
    model = model_functions[parsed_args["model"]]
    suffix = model_suffixes[parsed_args["model"]]
else
    error("--model argument must be one of $(keys(model_functions)); got $(parsed_args["model"])")
end

log_dN_default = make_log_dN(1.0, -1.0, 10.0, 0.7, 0.5)

samps = []
fnames = []
for fn in glob("GW*[0-9].h5", "/Users/wfarr/Research/o3a_posterior_samples/all_posterior_samples")
    h5open(fn, "r") do f
        push!(samps, read(f, "PublicationSamples/posterior_samples"))
        push!(fnames, fn)
    end
end
for fn in glob("*GW*_nocosmo.h5", "/Users/wfarr/Research/o3b_data/PE")
    h5open(fn, "r") do f
        push!(samps, read(f, "C01:Mixed/posterior_samples"))
        push!(fnames, fn)
    end
end

mask = map(samps) do ss
    median([x.mass_2_source for x in ss]) > mlow && median([x.mass_1_source for x in ss]) < mhigh
end
narrow_samps = samps[mask]
narrow_fnames = fnames[mask]
narrow_gwnames = [String(match(r"^.*(GW[0-9]+[_]*[0-9]*).*$", f)[1]) for f in narrow_fnames]

open(joinpath(@__DIR__, "..", "chains", "chain" * suffix * "_fnames.txt"), "w") do f
    for n in narrow_fnames
        write(f, n)
        write(f, "\n")
    end
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

m1s_sel, m2s_sel, zs_sel, pdraw, Ndraw = read_selection("/Users/wfarr/Research/o3b_data/O3-Sensitivity/endo3_bbhpop-LIGO-T2100113-v12.hdf5")
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
@info @sprintf("Minimum Neff_sel = %.1f, 4*nobs = %d", minimum(trace[:Neff_sel]), 4*sum(mask))
@info @sprintf("Minimum Neff_samps = %.1f", minimum([minimum(trace[n]) for n in namesingroup(trace, :Neff_samps)]))

idata = from_mcmcchains(trace, coords=Dict(:gwnames => narrow_gwnames), dims=(Neff_samps=(:gwnames,), m1s_popwt=(:gwnames,), m2s_popwt=(:gwnames,)), library="Turing")

to_netcdf(idata, joinpath(@__DIR__, "..", "chains", "chain" * suffix * ".nc"))