#!/bin/zsh

julia -t 4 sample.jl
julia -t 4 sample.jl --model two_broken_pl
julia -t 4 sample.jl --model broken_pl_plp
julia -t 4 sample.jl --model two_broken_pl_plp