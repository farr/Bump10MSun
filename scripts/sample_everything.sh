#!/usr/bin/env zsh
julia -t 4 sample.jl --model broken_pl --o4a true --mlow "$1"
julia -t 4 sample.jl --model broken_pl_gaussian --o4a true --mlow "$1"
julia -t 4 sample.jl --model broken_pl --o4a false --mlow "$1"
julia -t 4 sample.jl --model broken_pl_gaussian --o4a false --mlow "$1"

