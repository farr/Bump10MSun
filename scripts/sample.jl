using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Bump10MSun
using Glob
using HDF5
using Logging
using MCMCChainsStorage
using Printf
using StatsBase
using Turing

Nmcmc = 1000
Nchain = Threads.nthreads()
Npost = 128
Nsel = 4096

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

open(joinpath(@__DIR__, "..", "chains", "chains_fnames.txt"), "w") do f
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

model = broken_pl_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log.(pdraw), Ndraw)

if Nchain == 1
    trace = sample(model, NUTS(), Nmcmc)
else
    trace = sample(model, NUTS(), MCMCThreads(), Nmcmc, Nchain)
end
trace = with_logger(NullLogger()) do 
    append_generated_quantities(trace, generated_quantities(model, trace))
end
@info @sprintf("Minimum Neff_sel = %.1f, 4*nobs = %d", minimum(trace[:Neff_sel]), 4*sum(mask))
@info @sprintf("Minimum Neff_samps = %.1f", minimum([minimum(trace[n]) for n in namesingroup(trace, :Neff_samps)]))

h5open(joinpath(@__DIR__, "..", "chains", "chains.h5"), "w") do f
    write(f, trace)
end