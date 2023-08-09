#!/usr/bin/env zsh

julia -t 4 sample.jl --model broken_pl
julia -t 4 sample.jl --model broken_pl_gaussian

julia -t 4 sample.jl --model broken_pl --extra-events
julia -t 4 sample.jl --model broken_pl_gaussian --extra-events