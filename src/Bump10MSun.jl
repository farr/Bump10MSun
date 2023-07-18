module Bump10MSun

using Cosmology
using Distributions
using HDF5
using StatsBase
using StatsFuns
using Turing
using Unitful
using UnitfulAstro
using Trapz

include("model.jl")
export mlow, mhigh, zref
export MassFunction, BrokenPowerLaw, TwoBrokenPowerLaw, PowerLawGaussian
export PairingFunction, GaussianPairing, PowerLawPairing
export model_body
export make_log_dN, make_log_dN_tb
export broken_pl_model, two_broken_pl_model
export broken_pl_plp_model, two_broken_pl_plp_model
export power_law_plus_gaussian_model
export make_dNdm1, make_dNdm2, make_dNdq, make_dNdm, make_pairing_prob

include("weights.jl")
export li_prior_wt, md_sfr, md_sfr_zwt, read_selection
export resample_selection

include("utils.jl")
export with_seed, append_generated_quantities
export cumtrapz

end # module
