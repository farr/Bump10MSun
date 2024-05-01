module Bump10MSun

using ColorTypes
using Colors
using Cosmology
using Distributions
using GaussianKDEs
using Glob
using HDF5
using Interpolations
using LinearAlgebra
using Random
using StatsBase
using StatsFuns
using Turing
using Unitful
using UnitfulAstro
using Trapz

include("kde.jl")
export KDE, BoundedKDE
export logpdf_credible_levels

include("model.jl")
export mlow, mhigh, zref, mclow, mchigh
export m_lower_limit, m_upper_limit
export isselected, default_selection_fraction
export far_threshold
export MassFunction, BrokenPowerLaw, PowerLawGaussian
export PairingFunction, GaussianPairing, PowerLawPairing
export model_body
export make_log_dN, make_log_dN_tb
export broken_pl_model, two_broken_pl_model
export broken_pl_plp_model, two_broken_pl_plp_model
export power_law_plus_gaussian_model
export make_dNdm1, make_dNdm2, make_dNdq, make_dNdm, make_pairing_prob
export gaussian_model, two_gaussian_model

include("pputils.jl")
export median_plus_minus, median_hdi, result_macro
export suffix_map, label_map, var_name_map
export ns_var_names
export mf_label_map, mf_var_name_map
export pf_label_map, pf_var_name_map
export distribution_quantile
export bisect
export hdi_interval
export logscale_ticks
export categorical_palette

include("weights.jl")
export li_prior_wt, md_sfr, md_sfr_zwt, read_selection
export resample_selection

include("utils.jl")
export with_seed, append_generated_quantities
export cumtrapz
export load_pe_samples
export chirp_mass
export filter_selected

end # module
