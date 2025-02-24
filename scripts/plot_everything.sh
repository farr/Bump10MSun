#!/usr/bin/env zsh
julia plots_tables_macros.jl --o4a true --mlow "$1"
julia plots_tables_macros.jl --o4a false --mlow "$1"

