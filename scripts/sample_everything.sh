#!/usr/bin/env zsh

julia -t 4 sample.jl --model broken_pl_gaussian
julia -t 4 sample.jl --model broken_pl
julia -t 4 sample.jl --model gaussian

julia -t 4 sample.jl --model broken_pl_gaussian --selection_fraction 0.9
julia -t 4 sample.jl --model broken_pl --selection_fraction 0.9
julia -t 4 sample.jl --model gaussian --selection_fraction 0.9