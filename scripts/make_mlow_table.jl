using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Bump10MSun

Pkg.add(["DelimitedFiles", "Printf"])
using DelimitedFiles, Printf

data_dir = "../paper/figures/"
out_tex  = "../paper/table_2_content.tex"

name_bpl  = "m1pct_BrokenPL__including_230529"
name_bplg = "m1pct_BrokenPowerLaw+Gaussian__including_230529"

# ------------------------------------------------
# Keep only files like "<name>_<number>.txt"
# ------------------------------------------------
function files_with_numeric_suffix(data_dir::AbstractString, name::AbstractString)
    files = readdir(data_dir)
    out = String[]
    base = name * "_"
    for f in files
        if endswith(f, ".txt") && startswith(f, base)
            suffix = f[length(base)+1:end-4]
            if tryparse(Float64, suffix) !== nothing
                push!(out, f)
            end
        end
    end
    return sort(out)
end

# ------------------------------------------------
# Load grouped samples (same logic as your violin code)
# ------------------------------------------------
function load_grouped(data_dir::AbstractString, name::AbstractString)
    files = files_with_numeric_suffix(data_dir, name)

    samples = Dict{String, Vector{Float64}}()
    for file in files
        samples[file] = vec(readdlm(joinpath(data_dir, file)))
    end

    i = length(name * "_") + 1
    grouped = [samples[f] for f in files]
    mlows   = [parse(Float64, f[i:end-4]) for f in files]
    return mlows, grouped
end

# ------------------------------------------------
# Truncate (NOT round) to 2 decimals
# ------------------------------------------------
trunc2(x::Real) = trunc(x * 100) / 100
fmt2(x::Real) = @sprintf("%.2f", trunc2(x))

# ------------------------------------------------
# Extract HPD info
# Assumes: hpd(samples, alpha) -> [lo, mode, hi]
# alpha=0.10 -> 90% HPD
# ------------------------------------------------
function mode_pm_interval(samples::AbstractVector; alpha=0.10)
    h = hpd(samples, alpha)
    lo, mode, hi = h[1], h[2], h[3]
    plus  = hi - mode
    minus = mode - lo
    return lo, mode, hi, plus, minus
end

# ------------------------------------------------
# Load data
# ------------------------------------------------
mlows_bpl,  grouped_bpl  = load_grouped(data_dir, name_bpl)
mlows_bplg, grouped_bplg = load_grouped(data_dir, name_bplg)

bpl  = Dict(mlows_bpl[j]  => grouped_bpl[j]  for j in eachindex(mlows_bpl))
bplg = Dict(mlows_bplg[j] => grouped_bplg[j] for j in eachindex(mlows_bplg))

mlows = sort!(collect(intersect(keys(bpl), keys(bplg))))
isempty(mlows) && error("No overlapping mlow values between BPL and BPLG.")

# ------------------------------------------------
# Write LaTeX deluxetable with caption + subcolumn headers
# IMPORTANT: use raw strings so $ is not treated as interpolation
# ------------------------------------------------
caption = raw"\tablecaption{\textbf{Variation of $m_{1\%}$ as a function of $m_{\mathrm{low}}$. For each mass-model, we report the posterior mode and the 90\% highest density credible interval.}\label{tab:mlow-vs-m1}}"
head1   = raw"\colhead{$m_{\mathrm{low}} / M_\odot$} & \multicolumn{2}{c}{BPL} & \multicolumn{2}{c}{BPLG} \\"
head2   = raw" & \colhead{$m_{1\%} / M_\odot$ (90\%)} & \colhead{$m_{1\%} / M_\odot$ range (90\%)} & \colhead{$m_{1\%} / M_\odot$ (90\%)} & \colhead{$m_{1\%} / M_\odot$ range (90\%)} }"

open(out_tex, "w") do io
    println(io, raw"\begin{deluxetable}{lcccc}")
    println(io, caption)
    println(io, raw"\tablecolumns{5}")
    println(io, raw"\tablehead{")
    println(io, head1)
    println(io, head2)
    println(io, raw"\startdata")

    for m in mlows
        lo1, mode1, hi1, plus1, minus1 = mode_pm_interval(bpl[m])
        lo2, mode2, hi2, plus2, minus2 = mode_pm_interval(bplg[m])

        bpl_mode  = raw"$" * fmt2(mode1) * "^{+" * fmt2(plus1) * "}_{-" * fmt2(minus1) * "}" * raw"$"
        bpl_int   = "[" * fmt2(hi1) * ", " * fmt2(lo1) * "]"

        bplg_mode = raw"$" * fmt2(mode2) * "^{+" * fmt2(plus2) * "}_{-" * fmt2(minus2) * "}" * raw"$"
        bplg_int  = "[" * fmt2(hi2) * ", " * fmt2(lo2) * "]"

        println(io, fmt2(m), " & ", bpl_mode, " & ", bpl_int,
                    " & ", bplg_mode, " & ", bplg_int, raw" \\")
    end

    println(io, raw"\enddata")
    println(io, raw"\end{deluxetable}")
end

println("Wrote LaTeX table to: ", out_tex)
