using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ArviZ
using Bump10MSun
using CairoMakie
using CSV
using DataFrames
using DimensionalData
using Distributions
using LaTeXStrings
using NCDatasets
using Printf
using ProgressLogging
using Random
using StatsBase
using StatsFuns
using Trapz

# Set the theme---each of these "sections" and plots are wrapped in a begin ...
# end block so that you can run the whole thing within VSCode with a single
# Shift-Enter (execute the current block) command.  You can also run individual
# lines by selecting them with the mouse or keyboard and hitting Shift-Enter.
begin
    update_theme!(fontsize=12)
    update_theme!(pt_per_unit=0.75)
    colwidth_pt = 433.62
    colwidth = colwidth_pt / 0.75
    update_theme!(size=(colwidth,colwidth))
end

# Load the Data
begin
    sampso3a, fnameso3a, gwnameso3a = load_pe_samples("/Users/wfarr/Research/gwtc-2.1/", "*GW*_nocosmo.h5", "C01:Mixed/posterior_samples")
    sampso3b, fnameso3b, gwnameso3b = load_pe_samples("/Users/wfarr/Research/o3b_data/PE", "*GW*_nocosmo.h5", "C01:Mixed/posterior_samples")
    gwtc3_table = DataFrame(CSV.File("/Users/wfarr/Research/gwtc-3-tables/GWTC.csv"))

    samps = vcat(sampso3a, sampso3b)
    fnames = vcat(fnameso3a, fnameso3b)
    gwnames = vcat(gwnameso3a, gwnameso3b)

    samps, fnames, gwnames = filter_selected(samps, fnames, gwnames, gwtc3_table)
    @info "Analyzing $(length(gwnames)) events:"
    for n in gwnames
        @info "  $(n)"
    end
end

# Load the MCMC samples
traces = Dict(k => from_netcdf(joinpath(@__DIR__, "..", "chains", "chain" * suffix_map[k] * ".nc")) for k in keys(suffix_map))

# Compute various distributional quantities for the MCMC samples
begin
    ms = exp.(log(mlow):0.01:log(m_upper_limit))
    qs = collect(range(m_lower_limit/m_upper_limit, stop=1, length=128))
    ms3 = exp.(log(mlow):0.01:log(m_upper_limit))

    function make_distribution_map(dN_maker, xs)
        Dict(
            k => map([traces[k].posterior[v] for v in var_name_map[k[1:2]]]...) do R, args...
                dN = dN_maker(k[1], k[2], args...)
                R .* dN.(xs)
            end for k in keys(traces)
        )
    end

    dNdm1s = make_distribution_map(make_dNdm1, ms)
    dNdm2s = make_distribution_map(make_dNdm2, ms)
    dNdqs = make_distribution_map(make_dNdq, qs)

    pms = Dict(
        k => map([traces[k].posterior[v] for v in vcat(mf_var_name_map[k[1]], ns_var_names)]...) do R, args...
            p = make_dNdm(k[1], args...)
            pp = p.(ms)
            pp .= pp ./ trapz(ms, pp)
            pp
        end for k in keys(traces)
    )

    pms3 = Dict(
        k => map([traces[k].posterior[v] for v in vcat(mf_var_name_map[k[1]], ns_var_names)]...) do R, args...
            p = make_dNdm(k[1], args...)
            pp = p.(ms3)
            pp .= pp ./ trapz(ms3, pp)
            pp
        end for k in keys(traces)
    )

    m1pcts = Dict(
        k => [distribution_quantile(ms3, p, 0.01) for p in v]
        for (k,v) in pairs(pms3)
    )
end

# Figure 1
begin 
    Random.seed!(5766795581209736646)

    npt = 1024
    colors = categorical_palette(length(gwnames))
    f = Figure()
    a = Axis(f[1,1], xlabel=L"m_1 / M_\odot", ylabel=L"m_2 / M_\odot", limits=(0, m_upper_limit, 0, m_upper_limit))
    @progress for (i, (n, s)) in enumerate(zip(gwnames, samps))
        m1 = [x.mass_1_source for x in s]
        m2 = [x.mass_2_source for x in s]
        z = [x.redshift for x in s]
        wt = md_sfr_zwt.(z) ./ li_prior_wt.(m1, m2, z)
        inds = sample(1:length(z), Weights(wt), npt)
        m1 = m1[inds]
        m2 = m2[inds]

        pts = vcat(m1', m2')
        k = KDE(pts)
        inds = sample(1:size(pts, 2), 128)
        cpts = pts[:, inds]
        ls = logpdf_credible_levels(k, cpts, [0.1, 0.5])

        sig_m1 = std(m1)
        sig_m2 = std(m2)
            
        x = collect(range(minimum(m1)-sig_m1, stop=maximum(m1)+sig_m1, length=129))
        y = collect(range(minimum(m2)-sig_m2, stop=maximum(m2)+sig_m2, length=127))

        Z = [(xx > yy ? logpdf(k, [xx, yy]) : -Inf) for xx in x, yy in y]
        contour!(a, x, y, Z, levels=ls, color=colors[i], label=n)
    end

    axislegend(a, 
        [LineElement(color=colors[i]) for i in 1:length(gwnames)],
        gwnames,
        position=:lt,
        nbanks=3)

    band!(a, [mlow, m_upper_limit], [mlow, m_upper_limit], [m_upper_limit, m_upper_limit], color=(:grey, 0.5), label=nothing)
    band!(a, [0, mlow], [0,0], [m_upper_limit, m_upper_limit], color=(:grey, 0.5), label=nothing)

    m1s = collect(mhigh:0.01:m_upper_limit)
    m2s = [bisect(m2 -> chirp_mass(m1, m2)-mchigh, 0, m_upper_limit) for m1 in m1s]
    band!(a, m1s, m2s, m1s, color=(:grey, 0.5), label=nothing)

    save(joinpath(@__DIR__, "..", "paper", "figures", "m1-m2-contour.pdf"), f)
    f
end

# Figure 2
begin
    Random.seed!(9054612598149544883)

    colors = categorical_palette(2)
    xtickoptions = (xtickformat="{:.0f}", xticks=[3, 10, 20], xminorticks=range(mlow, mhigh, step=1), xminorticksvisible=true, xminorgridvisible=true)

    f = Figure()
    a1 = Axis(f[1,1]; xlabel=L"m_1 / M_\odot", ylabel=L"m_1 \mathrm{d} N / \mathrm{d} m_1 \mathrm{d} V \mathrm{d} t / \mathrm{Gpc}^{-3} \, \mathrm{yr}^{-1}", xscale=log10, yscale=log10, limits=(mlow, mhigh, 1e-1, 1e3), xtickoptions...)
    a2 = Axis(f[2,1]; xlabel=L"m_2 / M_\odot", ylabel=L"m_2 \mathrm{d} N / \mathrm{d} m_2 \mathrm{d} V \mathrm{d} t / \mathrm{Gpc}^{-3} \, \mathrm{yr}^{-1}", xscale=log10, yscale=log10, limits=(mlow, mhigh, 1e-2, 1e2), xtickoptions...)
    for (i, k) in enumerate([(BrokenPowerLaw(), PowerLawPairing()), (PowerLawGaussian(), PowerLawPairing())])
        lines!(a1, ms, ms .* mean(dNdm1s[k]), label=mf_label_map[k[1]], color=colors[i])
        for _ in 1:100
            lines!(a1, ms, ms .* sample(dNdm1s[k]), color=(colors[i], 0.1), label=nothing)
        end

        lines!(a2, ms, ms .* mean(dNdm2s[k]), color=colors[i], label=nothing)
        for _ in 1:100
            lines!(a2, ms, ms .* sample(dNdm2s[k]), color=(colors[i], 0.1), label=nothing)
        end
    end

    axislegend(a1, position=:rb)

    save(joinpath(@__DIR__, "..", "paper", "figures", "dNdm_traces.pdf"), f)
    f
end

# Figure 3
begin 
    Random.seed!(3495738867787726282)

    yl, yh = 1e-2, 10
    colors = categorical_palette(2)

    f = Figure()
    a = Axis(f[1,1], xlabel=L"m / M_\odot", ylabel=L"m p(m)", limits=(mlow, mhigh, yl, yh), xscale=log10, yscale=log10, xticks=[3, 10, 20], xtickformat="{:.0f}", xminorticks=range(3, 20, step=1), xminorgridvisible=true, xminorticksvisible=true)

    for (i, k) in enumerate([(BrokenPowerLaw(), PowerLawPairing()), (PowerLawGaussian(), PowerLawPairing())])
        pm = pms[k]

        lines!(a, ms, ms .* mean(pm), label=mf_label_map[k[1]], color=colors[i])
        for _ in 1:100
            lines!(a, ms, ms .* sample(pm), label=nothing, color=(colors[i], 0.1))
        end
    end

    axislegend(a, position=:lt)

    save(joinpath(@__DIR__, "..", "paper", "figures", "pm_traces.pdf"), f)
    f
end

# Figure 4
begin 
    xs = 3:0.01:6

    f = Figure()
    a = Axis(f[1,1], xlabel=L"m_{1%} / M_\odot", limits=(minimum(xs), maximum(xs), 0, nothing), xminorgridvisible=true, yminorgridvisible=true, palette=(color=categorical_palette(2),))

    for (i, k) in enumerate([(BrokenPowerLaw(), PowerLawPairing()), (PowerLawGaussian(), PowerLawPairing())])
        sf = default_selection_fraction
        kde = BoundedKDE(vec(m1pcts[k]), lower=3)
        lines!(a, xs, pdf.((kde,), xs), label=mf_label_map[k[1]])
    end

    axislegend(a, position=:rt)

    save(joinpath(@__DIR__, "..", "paper", "figures", "m1pct.pdf"), f)
    f
end

# Macros
begin
    function write_macro(f, s)
        @info "Writing macro: $(s)"
        write(f, s)
    end

    open(joinpath(@__DIR__, "..", "paper", "result_macros.tex"), "w") do f
        write_macro(f, result_macro(raw"\dNlogmpeak", raw"\mathrm{Gpc}^{-3} \, \mathrm{yr}^{-1}", traces[PowerLawGaussian(), PowerLawPairing()].posterior.R, digits=0))
        write_macro(f, result_macro(raw"\monepctplgplp", raw"M_\odot", m1pcts[PowerLawGaussian(), PowerLawPairing()], digits=2))
        write_macro(f, result_macro(raw"\mpeakplgplp", raw"M_\odot", traces[PowerLawGaussian(), PowerLawPairing()].posterior.mu, digits=2))
        write_macro(f, result_macro(raw"\alphatwoplgplp", traces[PowerLawGaussian(), PowerLawPairing()].posterior.a2, digits=1))
        write_macro(f, "\\newcommand{\\mlow}{$(mlow)}\n\\newcommand{\\mlowunits}{\\ensuremath{\\mlow \\, M_\\odot}}")
        write_macro(f, "\\newcommand{\\mclow}{$(mclow)}\n\\newcommand{\\mclowunits}{\\ensuremath{\\mclow \\, M_\\odot}}")
        write_macro(f, "\\newcommand{\\mchigh}{$(mchigh)}\n\\newcommand{\\mchighunits}{\\ensuremath{\\mchigh \\, M_\\odot}}")
        write_macro(f, "\\newcommand{\\mhigh}{$(mhigh)}\n\\newcommand{\\mhighunits}{\\ensuremath{\\mhigh \\, M_\\odot}}")
        write_macro(f, "\\newcommand{\\nevts}{$(length(gwnames))}")
    end
end

# Table 1
begin
    open(joinpath(@__DIR__, "..", "paper", "table1_content.tex"), "w") do f
        write(f, "\\begin{deluxetable}{llll}\n\\tablecolumns{3}\n\\tablecaption{\\label{tab:monepct} \$m_{1\\%}\$ for our various models and using different selection functions.}\n")
        write(f, "\\tablehead{\\colhead{Mass Function Model} & \\colhead{\$m_{1\\%} / M_\\odot\$ (1-\$\\sigma\$, 68\\%)} & \\colhead{\$m_{1\\%} / M_\\odot\$ range (2-\$\\sigma\$, 95\\%)}}\n")
        write(f, "\\startdata\n")
        for mf in [BrokenPowerLaw(), PowerLawGaussian()]
            k = (mf, PowerLawPairing())
            m1p = m1pcts[k]
    
            m, h, l = median_plus_minus(vec(m1p), 0.68)
            ll, hh = quantile.((vec(m1p),), [0.025, 0.975])
    
            write(f, "\\texttt{" * mf_label_map[mf] * "}")
            write(f, @sprintf("& \$%.2f^{+%.2f}_{-%.2f}\$", m, h-m, m-l))
            write(f, @sprintf("& \$\\left[%.1f, %.1f \\right]\$", ll, hh))
            write(f, "\\\\ \n")
        end
        write(f, "\\enddata\n\\end{deluxetable}\n")
    end
end