using Pkg
Pkg.gc()
Pkg.add(["Colors", "Plots", "StatsPlots", "LaTeXStrings"])
using DelimitedFiles, StatsPlots, LaTeXStrings

# Directory containing the files
data_dir = "../paper/figures/"  # Change this to your actual directory

# Get all BPL files in the directory
name = "m1pct_BrokenPL__including_230529"
files = filter(f -> occursin(name, f) && endswith(f, ".txt"), readdir(data_dir))

# Read data from each file
samples = Dict()
for file in files
    filepath = joinpath(data_dir, file)
    samples[file] = vec(readdlm(filepath))  # Read file and flatten into a vector
end


# Convert dictionary to a format suitable for StatsPlots
i= length(name * "_") + 1
grouped_data_BPL = [samples[f] for f in files]
flat_data = vcat(grouped_data_BPL...)
println([f[i:length(f)-4] for f in files])
labels_BPL = [[parse(Float64, f[i:length(f)-4])].*ones(length(grouped_data_BPL[j])) for (j,f) in enumerate(files)] # X-axis labels (filenames)

# Get all text files in the directory
name = "m1pct_BrokenPowerLaw+Gaussian__including_230529"
files = filter(f -> occursin(name, f) && endswith(f, ".txt"), readdir(data_dir))

# Read data from each file
samples = Dict()
for file in files
    filepath = joinpath(data_dir, file)
    samples[file] = vec(readdlm(filepath))  # Read file and flatten into a vector
end


# Convert dictionary to a format suitable for StatsPlots
i= length(name * "_") + 1
grouped_data_BPLG = [samples[f] for f in files]
flat_data = vcat(grouped_data_BPLG...)
println([f[i:length(f)-4] for f in files])
labels_BPLG = [[parse(Float64, f[i:length(f)-4])].*ones(length(grouped_data_BPLG[j])) for (j,f) in enumerate(files)] # X-axis labels (filenames)

# Create a violin plot
# Define LCHuv colors
lch_colors = [LCHuv(50.0, 90.0, 90.0), LCHuv(50.0, 90.0, 270.0)]

# Convert to RGB format
rgb_colors = [convert(RGB, c) for c in lch_colors]
transparent_colors = [RGBA(c.r, c.g, c.b, 0.5) for c in rgb_colors]  # Adjust alpha as needed

colors = [vcat([[transparent_colors[1] for k in 1:length(grouped_data_BPL[j])] for j in 1:length(files)]...), vcat([[transparent_colors[2] for k in 1:length(grouped_data_BPLG[j])] for j in 1:length(files)]...) ]

violin(labels_BPL, grouped_data_BPL, ylabel=L"m_{1\%}/M_{\odot}", xlabel=L"m_\mathrm{low}/M_{\odot}", title="", color = colors[1], framestyle = :box, label="")#, bandwidth = 1.0)#, label = ["BrokenPL", "", "", "", ""])#, ylabelfontsize=21, ylim=(2.0, 5.0))
violin!(labels_BPLG, grouped_data_BPLG, ylabel=L"m_{1\%}/M_{\odot}", xlabel=L"m_\mathrm{low}/M_{\odot}", title="", color = colors[2], framestyle = :box, label="")#, bandwidth=2.0)#, label = ["BrokenPL+G", "", "", "", ""])#, ylabelfontsize=21, ylim=(2.0, 5.0))
scatter!([NaN], [NaN], color=colors[1], label="BPL", lw=3, markershape=:square, markersize=8)
scatter!([NaN], [NaN], color=colors[2], label="BPLG", lw=3, markershape=:square, markersize=8)
plot!(legend=:topright)
# Show plot
savefig("../paper/figures/violin_plot.pdf")  # Save if needed
