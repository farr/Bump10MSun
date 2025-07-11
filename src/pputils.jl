using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Interpolations
using Statistics
using KernelDensity

"""
    hdi_interval(xs::AbstractVector, q)

Return `(l,h)`, the shortest interval that contains a fraction `q` of the values
in `xs`.
"""

function hdi_interval(xs::AbstractVector, q)
    xs = sort(xs)
    n = length(xs)
    nn = round(Int, n*q)

    dx = xs[nn:end] .- xs[1:end-nn+1]
    i = argmin(dx)

    (xs[i], xs[i+nn-1])
end

"""
    median_plus_minus(x, q=0.68)

Return the median and `(q+1)/2` and `(1-q)/2` quantiles of the collection `x`.
"""
function median_plus_minus(x, q=0.68)
    @assert (zero(q) < q && q < one(q)) "quantile `q` must be between zero and one"
    x = vec(x)
    m = median(x)
    h = quantile(x, (q+1)/2)
    l = quantile(x, (1-q)/2)
    (m, h, l)
end

"""
    median_hdi(x, q=0.68)

Return median and highest density `q` interval of the collection `x`.
"""
function median_hdi(x, q=0.68)
    @assert (zero(q) < q && q < one(q))
    x = vec(x)
    m = median(x)
    l, h = hdi_interval(x, q)

    return (m, h, l)
end

function hpd(samples, alpha)
    x = sort(copy(samples))  # Sort a copy of the samples
    n = length(x)
    cred_mass = 1.0 - alpha

    interval_idx_inc = floor(Int, cred_mass * n)
    n_intervals = n - interval_idx_inc
    interval_width = x[interval_idx_inc+1:end] .- x[1:n_intervals]

    if length(interval_width) == 0
        throw(ArgumentError("Too few elements for interval calculation"))
    end

    min_idx = argmin(interval_width)
    hdi_min = x[min_idx]
    hdi_max = x[min_idx + interval_idx_inc]

    kde_result = kde(samples)
    # Find the mode by locating the maximum density value
    mode_index = argmax(kde_result.density)
    mode_approx = kde_result.x[mode_index]

    return [hdi_min, mode_approx, hdi_max]
end



raw"""
    result_macro(mname, samples; digits=1, q=0.68)
    result_macro(mname, units, samples; digits=1, q=0.68)

Return a string with a LaTeX macro defining a credible `q` range based on `samples`.

Round macro values using `round` to the given `digits`.  If given `units` produce an additional 
LaTeX macro with `$(mname)units` that gives the values and includes the unit commands given in `units`.

# Examples
```jldoctest
julia> result_macro(raw"\height", raw"\mathrm{m}", range(1.9, stop=2.1, length=201), digits=3)
"\\newcommand{\\height}{\\ensuremath{2.0^{+0.068}_{-0.068}}}\n\\newcommand{\\height_units}{\\ensuremath{\\height \\, \\mathrm{m}}}"
```
"""
function result_macro(mname, samples; digits=1, q=0.68)
    m,h,l = median_plus_minus(samples, q)
    if digits <= 0
        mr = Int(round(m, digits=digits))
        hr = Int(round(h-m, digits=digits))
        lr = Int(round(m-l, digits=digits))
    else
        mr = round(m, digits=digits)
        hr = round(h-m, digits=digits)
        lr = round(m-l, digits=digits)
    end
    "\\newcommand{$(mname)}{\\ensuremath{$(mr)^{+$(hr)}_{-$(lr)}}}\n"
end
function result_macro(mname, units, samples; digits=1, q=0.68)
    result_macro(mname, samples; digits=digits, q=q) * "\\newcommand{$(mname)units}{\\ensuremath{$(mname) \\, $(units)}}\n"
end

"""Map from model to suffix for chain files"""
suffix_map = Dict(
    (BrokenPowerLaw(), PowerLawPairing()) => "_bpl",
    (PowerLawGaussian(), PowerLawPairing()) => "_bplg"
)

"""Map from model to variable names for mass function variables."""
mf_var_name_map = Dict(
    BrokenPowerLaw() => [:R, :a1, :a2, :mb],
    PowerLawGaussian() => [:R, :a1, :a2, :mu, :sigma, :fg]
)

"""Map from model to variable names for pairing function variables."""
pf_var_name_map = Dict(
    GaussianPairing() => [:mu_q, :sigma_q],
    PowerLawPairing() => [:beta]
)

ns_var_names = [:r_ns, :mu_ns, :sigma_ns]

"""Map from model to variable names for all variables."""
var_name_map = Dict((km, kp) => vcat(mv, pv, ns_var_names) for (km, mv) in pairs(mf_var_name_map) for (kp, pv) in pairs(pf_var_name_map))

"""Map from model to plot labels for mass function models."""
mf_label_map = Dict(
    BrokenPowerLaw() => "Broken PL",
    PowerLawGaussian() => "Broken Power Law + Gaussian"
)
"""Map from model to plot labels for pairing function models."""
pf_label_map = Dict(
    GaussianPairing() => "Gaussian Pair",
    PowerLawPairing() => "Power Law Pair"
)

"""Map from model to plot labels for combined models."""
label_map = Dict(
    k => mf_label_map[k[1]] * ", " * pf_label_map[k[2]] for k in keys(suffix_map)
)

"""
    distribution_quantile(f, xs, q)
    distribution_quantile(xs, ys, q)

Return an estimate of the `x` that corresponds to the `q` quantile of the
distribution `f(x)` or `ys`.
"""
function distribution_quantile(f::Function, xs::AbstractVector, q::Real)
    fm = f.(xs)
    distribution_quantile(xs, fm, q)
end
function distribution_quantile(xs::AbstractVector, ys::AbstractVector, q::Real)
    cfm = cumtrapz(xs, ys)
    cfm .= cfm ./ cfm[end]
    LinearInterpolation(cfm, xs)(q)
end

"""
    bisect(f, xl, xh)

Return a solution to `f(x) = 0` between `xl` and `xh` using the bisection method.
"""
function bisect(f, xl, xh)
    bisect_loop(f, f(xl), f(xh), xl, xh)
end
function bisect_loop(f, fl, fh, xl, xh)
    if xh - xl < 1e-8
        (xh+xl)/2
    else
        xm = (xh+xl)/2
        fm = f(xm)

        if fl*fm < 0
            bisect_loop(f, fl, fm, xl, xm)
        else
            bisect_loop(f, fm, fh, xm, xh)
        end
    end
end

"""
    bounded_kde_pdf(xs::AbstractArray; low=nothing, high=nothing)

Return a PDF function for a KDE of `xs` on a domain that is bounded by `low` and
`high`.
"""
function bounded_kde_pdf(xs::AbstractArray; low=nothing, high=nothing)
    k = KDE(xs)
    function mypdf(x)
        p = pdf(k, [x])
        if low !== nothing
            p += pdf(k, [2*low-x])
        end
        if high !== nothing
            p += pdf(k, [2*high-x])
        end
        p
    end
    mypdf
end

"""
    categorical_palette(n; luminence=50, chroma=90)

Return a palette of `n` distinguishable colors.
"""
function categorical_palette(n; luminence=50, chroma=90)
    [LCHuv(luminence, chroma, 360*(i-1/2)/n) for i in 1:n]
end
