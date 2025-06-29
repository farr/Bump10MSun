using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.add(PackageSpec(name="StatsPlots", version="0.14.33"))

using DelimitedFiles, DataFrames, StatsPlots, LaTeXStrings, Plots

# Load from text files (assumes one column of numbers per file)
m1pct = readdlm(joinpath(@__DIR__, "..", "paper", "figures","m1pct_BPLG__including_230529_3.0.txt"))[:]
mu = readdlm(joinpath(@__DIR__, "..", "paper", "figures","mu_BPLG__including_230529_3.0.txt"))[:]
sigma = readdlm(joinpath(@__DIR__, "..", "paper", "figures","sigma_BPLG__including_230529_3.0.txt"))[:]

# Assemble into DataFrame with LaTeX-style labels
df = DataFrame(L"m_{1\%}" => m1pct, L"m_{peak}" => mu, L"\sigma" => sigma)

# Generate corner plot
p = @df df StatsPlots.pairplot([L"m_{1\%}", L"m_{peak}", L"\sigma"], label = "", grid = false)
savefig(p, joinpath(@__DIR__, "..", "paper", "figures", "correlations_3.pdf"))
