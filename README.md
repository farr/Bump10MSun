# Fitting Broken Power Laws to the Low-Mass LVK Bump

[![Paper PDF](https://img.shields.io/badge/Paper-PDF-blue)](https://github.com/farr/Bump10MSun/blob/main-pdf/paper/Bump10MSun.pdf)

What it says on the tin.  We only fit the O3 confident events (O3a and O3b)
because these are the only events for which there are actual injections used to
estimate sensitivity.  (O1 and O2 used a semi-analytic approximation for the
detector sensitivity.)

If you want to reproduce the results here you will need to download the O3a
posterior samples [here](https://dcc.ligo.org/LIGO-P2000223-v7/public), the O3b
posterior samples [here](https://zenodo.org/record/5546663#.Yf9Uge7ML0o), and
the search sensitivity data files for O3
[here](https://zenodo.org/record/5546676#.Yf9UMe7ML0o) and modify the paths in
`scripts/sample.jl` accordingly.

## Installation

Install Julia from [julialang.org](https://julialang.org/downloads/).

Clone this repository.

Instantiate the necessary Julia environment: 

1. Start Julia.

2. Instantiate the Julia environment.  For now you will also have to manually
   (re)add the [PopModels.jl](https://github.com/farr/PopModels.jl) package, as
   it is not yet listed in the General Julia registry.  Push `]` to enter the
   package manager mode, activate the current directory with `activate .`, and
   issue the command `add https://github.com/farr/PopModels.jl`.  Now that the
   PopModels.jl package is registered, you can instantiate the environment by
   issuing the `instantiate` command.

3. You should now have a working `Bump10MSun` environment, and can:
    - Re-run the MCMC sampling of the various population models using the
      `scripts/sample.jl` script (artisinal) or resample all the models using
      `scripts/sample_everything.sh`.
    - Examine and adapt the population models in `src/model.jl` if you would
      prefer to change the "common" mass function or the pairing function.
    - Regenerate plots numbers, etc that appear in the paper draft / sketch
      using the notebook at `notebooks/PaperPlotsTablesMacros.ipynb`.
    - Edit and/or recompile the paper LaTeX in `paper` using `latexmk` or
      similar.  (Note that a PDF of the paper is compiled automatically at
      GitHub via CI actions if you push a branch.)