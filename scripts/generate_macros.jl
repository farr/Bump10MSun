using ArviZ
using Bump10MSun

ms = exp.(range(log(3), stop=log(20), length=128))

traces = Dict(k => from_netcdf(joinpath(@__DIR__, "..", "chains", "chain" * suffix_map[k] * ".nc")) for k in keys(suffix_map))

pm_map = Dict(
    k => map([traces[k].posterior[v] for v in mf_var_name_map[k[1]]]...) do args...
        make_dNdm(k[1], args[2:end]...)
    end for k in keys(traces)
)

function write_macro(f, s)
    @info "Writing macro: $(s)"
    write(f, s)
end

open(joinpath(@__DIR__, "..", "paper", "result_macros.tex"), "w") do f
    write_macro(f, result_macro(raw"\dNlogmpeak", raw"\mathrm{Gpc}^{-3} \, \mathrm{yr}^{-1}", traces[PowerLawGaussian(), PowerLawPairing()].posterior.R, digits=-1))
    write_macro(f, result_macro(raw"\monepctplgplp", raw"M_\odot", distribution_quantile.(pm_map[PowerLawGaussian(), PowerLawPairing()], (ms,), 0.01), digits=1))
    write_macro(f, result_macro(raw"\mpeakplgplp", raw"M_\odot", traces[PowerLawGaussian(), PowerLawPairing()].posterior.mu, digits=2))
    write_macro(f, result_macro(raw"\alphatwoplgplp", traces[PowerLawGaussian(), PowerLawPairing()].posterior.a2, digits=1))
    write_macro(f, "\\newcommand{\\mclow}{$(mclow)}\n\\newcommand{\\mclowunits}{\\ensuremath{\\mclow \\, M_\\odot}}")
    write_macro(f, "\\newcommand{\\mchigh}{$(mchigh)}\n\\newcommand{\\mchighunits}{\\ensuremath{\\mchigh \\, M_\\odot}}")
end