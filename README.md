# Fitting Broken Power Laws to the Low-Mass LVK Bump

[![Paper PDF](https://img.shields.io/badge/Paper-PDF)](https://github.com/farr/Bump10MSun/blob/main-pdf/paper/Bump10MSun.pdf)

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
