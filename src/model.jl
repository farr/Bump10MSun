"""The minimum mass for plotting."""
const mlow = 3.0

"""The lower limit of the mass function used in integrals."""
const m_lower_limit = 1.0

"""The upper limit of the mass function used in integrals."""
const m_upper_limit = 50.0

"""The maximum mass for plotting."""
const mhigh = 20.0

"""The minimum chirp mass considered in this study."""
const mclow = 2.612

"""The maximum chirp mass considered in this study."""
const mchigh = 17.41

"""The reference redshift, at which rates are evaluated."""
const zref = 0.0

"""If we need a mass reference, we use this."""
const mref = 10.0

"""The fraction of likelihood-weighted PE samples that must pass the selection
criterion to be included in this study."""
const default_selection_fraction = 0.5

"""Below this mass, use NS mass function; above, BH."""
const mmax_ns = 3.0

square(x) = x*x

# We use these types as tags to ensure that we evaluate the right functions,
# since they cannot be distinguished by argument number.
abstract type MassFunction end
abstract type PairingFunction end

struct BrokenPowerLaw <: MassFunction end
struct PowerLawGaussian <: MassFunction end

struct GaussianPairing <: PairingFunction end
struct PowerLawPairing <: PairingFunction end

function isselected(m1, m2)
  mlow <= m1 && chirp_mass(m1, m2) <= mchigh
end

# It is convenient that all of these different models are distinguished by their number of arguments.
function _fm(::BrokenPowerLaw, m, a1, a2, mb, fns, mu_ns, sigma_ns)
  x = m/mb
  if m < mmax_ns
    log(fns) - 0.5*square((m-mu_ns)/sigma_ns)
  elseif m < mb
      a1*log(x)
  else
      a2*log(x)
  end
end

function _fm(::PowerLawGaussian, m, a1, a2, mu, sigma, fg, fns, mu_ns, sigma_ns)
  if m < mmax_ns
    log(fns) - 0.5*square((m-mu_ns)/sigma_ns)
  else
    log_fg = log(fg)
    log_fpl = log1p(-fg)
    if m < mu
      log_plg = logaddexp(log_fpl + a1*log(m/mu), log_fg - 0.5*square((m-mu)/sigma))
    else
      log_plg = logaddexp(log_fpl + a2*log(m/mu), log_fg - 0.5*square((m-mu)/sigma))
    end
  end
end

_pf(::GaussianPairing, q, mu, sigma) = -0.5*square((q-mu)/sigma) + 0.5*square(1-mu)/sigma
_pf(::PowerLawPairing, q, beta) = beta*(log1p(q) - log(2))

"""
  make_log_dN(BrokenPowerLaw(), GaussianPairing(), a1, a2, mb, mu_q, sigma_q, fns, mu_ns, sigma_ns)

Produces a function that computes `(m1, m2) -> log_dNdm1dm2` for the
single-break power law population model.

The population model is 

```math
m_1 m_2 \\frac{\\mathrm{d}N}{\\mathrm{d} m_1 \\mathrm{d} m_2 \\mathrm{d} V \\mathrm{d} t} = R f_m\\left( m_1 \\right) f_m\\left( m_2 \\right) \\exp\\left( -\\frac{\\left(q - \\mu_q\\right)^2}{2 \\sigma_q^2} + \\frac{\\left(1-\\mu_q\\right)^2}{2 \\sigma_q^2} \\right)^2
```
where
```math
f_m\\left( m \\right) = \\begin{cases}
m^{\\alpha_1} & m_\\mathrm{min} < m < m_b \\
m^{\\alpha_2} & m_b \\leq m < m_\\mathrm{max}
\\end{cases}
```
is the "common" black hole mass function, and a Gaussian-shaped pairing function
in ``q`` is used with parameters ``mu_q`` and ``sigma_q``.  The function
returned is normalized so that ``R`` measures the volumetric merger rate per
natural log mass squared at `m1 = m2 = m_b`.

TODO: document NS parameters.
"""
function make_log_dN(mf::BrokenPowerLaw, pf::GaussianPairing, a1, a2, mb, mu_q, sigma_q, r_ns, mu_ns, sigma_ns)
function log_dN(m1, m2)
    if isselected(m1, m2)
      q = m2/m1
      _fm(mf, m1, a1, a2, mb, r_ns, mu_ns, sigma_ns) + _fm(mf, m2, a1, a2, mb, r_ns, mu_ns, sigma_ns) - log(m1) - log(m2) + _pf(pf, q, mu_q, sigma_q)
    else
      -Inf
    end
end
log_dN
end

"""
  make_log_dN(BrokenPowerLaw(), PowerLawPairing(), a1, a2, mb, beta, fns, mu_ns, sigma_ns)

Produces a function that computes `(m1, m2) -> log_dNdm1dm2` for the
single-break power law population model.

The population model is 

```math
m_1 m_2 \\frac{\\mathrm{d}N}{\\mathrm{d} m_1 \\mathrm{d} m_2 \\mathrm{d} V \\mathrm{d} t} = R f_m\\left( m_1 \\right) f_m\\left( m_2 \\right) \\left(\\frac{1 + q}{2}\\right)^\\beta
```
where
```math
f_m\\left( m \\right) = \\begin{cases}
m^{\\alpha_1} & m_\\mathrm{min} < m < m_b \\
m^{\\alpha_2} & m_b \\leq m < m_\\mathrm{max}
\\end{cases}
```
is the "common" black hole mass function, and a power-law pairing function in
``1+q`` is used with parameter ``beta``.  The function returned is normalized so
that ``R`` measures the volumetric merger rate per natural log mass squared at
`m1 = m2 = m_b`.

TODO: document NS parameters.
"""
function make_log_dN(mf::BrokenPowerLaw, pf::PowerLawPairing, a1, a2, mb, beta, r_ns, mu_ns, sigma_ns)
  function log_dN(m1, m2)
    if isselected(m1, m2)
      q = m2/m1
      _fm(mf, m1, a1, a2, mb, r_ns, mu_ns, sigma_ns) + _fm(mf, m2, a1, a2, mb, r_ns, mu_ns, sigma_ns) - log(m1) - log(m2) + _pf(pf, q, beta)
    else
      -Inf
    end
  end
  log_dN
end

"""
  make_log_dN(PowerLawGaussian(), PowerLawPairing(), a1, a2, mu, sigma, fg,
  beta)

Produces a function that computes `(m1, m2) -> log_dNdm1dm2` for the power-law
plus Gaussian population model.

The population model is 
```math
m_1 m_2 \\frac{\\mathrm{d}N}{\\mathrm{d} m_1 \\mathrm{d} m_2 \\mathrm{d} V \\mathrm{d} t} = R f_m\\left( m_1 \\right) f_m\\left( m_2 \\right) \\left( \\frac{1+q}{2} \\right)^\\beta
```
where
```math
f_m \\left( m \\right) = f_g \\exp\\left(-\\frac{1}{2} \\left(\\frac{m-\\mu}{\\sigma} \\right)^2 \\right) + \\left( 1 - f_g \\right) \\left( \\frac{m}{\\mu} \\right)^{\\alpha_1, \\alpha_2}
```
is a sum of broken power-law "continuium" and a Gaussian "peak."  ``\\alpha_1``
is used when ``m < \\mu`` and ``\\alpha_2`` when ``m > \\mu``.  As written, `R`
is the merger rate per log mass squared at `m1 = m2 = mu`.
"""
function make_log_dN(mf::PowerLawGaussian, pf::PowerLawPairing, a1, a2, mu, sigma, fg, beta, r_ns, mu_ns, sigma_ns)
function log_dN(m1, m2)
  if isselected(m1, m2)
    q = m2/m1
    _fm(mf, m1, a1, a2, mu, sigma, fg, r_ns, mu_ns, sigma_ns) + _fm(mf, m2, a1, a2, mu, sigma, fg, r_ns, mu_ns, sigma_ns) - log(m1) - log(m2) + _pf(pf, q, beta)
  else
    -Inf
  end
end
log_dN
end

"""
  make_dNdm1(mf_type, pair_type, args...)

Returns a function that computes the marginal density (unnormalized) in `m1`.
Mass function and pairing function types and arguments are passed to ``make_log_dN``.  
"""
function make_dNdm1(mf_type, pair_type, args...)
  log_dN = make_log_dN(mf_type, pair_type, args...)
  make_dNdm1(log_dN)
end
function make_dNdm1(log_dN)
  m1 -> begin
      if m1 <= m_lower_limit
          zero(m1)
      elseif m1 > m_upper_limit
          zero(m1)
      else
          n = Int(round(100*(log(m1)-log(m_lower_limit)))) + 2
          ms = exp.(range(log(m_lower_limit), stop=log(m1), length=n))
          dN = exp.(map(m2->log_dN(m1,m2), ms))
          trapz(ms, dN)
      end
  end
end

"""
  make_dNdm2(mf_type, pair_type, args...)

Returns a function that computes the marginal density (unnormalized) in `m2`;
mass function and pairing function types and ``args...`` are passed to
``make_log_dN``.
"""
function make_dNdm2(mf_type, pair_type, args...)
  log_dN = make_log_dN(mf_type, pair_type, args...)
  make_dNdm2(log_dN)
end
function make_dNdm2(log_dN)
  m2 -> begin
      if m2 < m_lower_limit
          zero(m2)
      elseif m2 > m_upper_limit
          zero(m2)
      else
          n = Int(round(100*(log(m_upper_limit) - log(m2)))) + 2
          ms = exp.(range(log(m2), stop=log(m_upper_limit), length=n))
          dN = exp.(map(m1->log_dN(m1, m2), ms))
          trapz(ms, dN)
      end
  end
end

"""
Returns a function that computes the marginal (unnormalized) density in `q`.
Mass function type, pairing function type, and arguments are passed to
`make_log_dN`.
"""
function make_dNdq(mf_type, pf_type, args...)
  log_dN = make_log_dN(mf_type, pf_type, args...)
  make_dNdq(log_dN)
end
function make_dNdq(log_dN)
  q -> begin
      if q < m_lower_limit/m_upper_limit
          zero(q)
      elseif q > 1
          zero(q)
      else
          ml = m_lower_limit/q
          mh = m_upper_limit
          n = Int(round(100*(log(mh)-log(ml)))) + 2
          ms = exp.(range(log(ml), stop=log(mh), length=n))
          dN = exp.(map(m1->log_dN(m1, q*m1) + log(m1), ms))
          trapz(ms, dN)
      end
  end
end

"""
Returns a function representing the "common" black hole mass function, ``f_m``
above; with two slopes and one break, the single broken power law; with three
slopes and two breaks, the double-broken power law.
"""
function make_dNdm(mf_type, args...)
m -> begin
  exp(_fm(mf_type, m, args...))/m
end
end

"""
Returns the "pairing probability", the function ``g(q)``.
"""
function make_pairing_prob(pf_type, args...)
  q -> exp(_pf(pf_type, q, args...))
end

"""
  model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)

Implements the mass model fitting, given a log density function `log_dN` and
observed masses (drawn from a prior with weight `log_wts`) and detected
injection masses with (normalized) draw distribution `log_pdraw`.  

# Arguments

- `log_dN``: log density function.  Called as `log_dN(m1, m2)`.
- `m1s``: observed primary masses, nested iterables of lengeth `n_observations`,
then `n_samples` within each observation.
- `m2s``: observed secondary masses.
- `log_wts`: log of the prior weight assigned to each observed mass pair.
- `m1s_sel``: primary masses from detected injections, iterable of length
`n_detected`.
- `m2s_sel`: secondary masses from detected injections.
- `log_pdraw`: log of the normalized draw distribution for the injections.
- `Ndraw`: the total number of injections drawn to produce the detections in
`m1s_sel`, `m2s_sel`.

# Returns

`(log_likelihood_sum, log_normalization_sum, generated_quantities)`

After calling `model_body(...)`, you should add the log-likelihood and
log-normalization terms to the density via

  Turing.@addlogprob! log_likelihood_sum
  Turing.@addlogprob! log_normalization_sum

The generated quantities include 

- `R`: The merger rate density scaling parameter (overall merger rate density is
`R*exp(log_dN(...))`).
- `Neff_sel`: The number of effective samples in the importance weighted
integral of the selection function (the normalization); see [Farr
(2019)](https://ui.adsabs.harvard.edu/abs/2019RNAAS...3...66F/abstract) for
definition.  This should be at absolute minimum ``4 N_\\mathrm{obs}`` (i.e. 4
times the catalog size), preferably much larger.
- `Neff_samps`: An array giving the number of effective samples in the
importance-weighted likelihood integral for each observation; these should all
be ``\\gg 1`` (3-4 is good, 10 is better) or else the likelihood integral is
not converged.
- `m1s_popwt`: An array giving a population-weighted draw from the primary mass
samples for each event.
- `m2s_popwt`: Same, for secondary masses.
- `m1_draw`: A draw from the primary masses of detected injections weighted by
the population (these draws can be used to predict the *observed* population
implied by the model).
- `m2_draw`: Same, but secondary masses.

You should return the generated quantities from the turing model, so the full
use of this function looks like 

  @model function my_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw, ...)
    # Set up parameter priors, etc.
    my_parameter_1 ~ Normal(...)
    # and so on

    # Compute any derived quantities you need.
    derived_quantity = my_parameter_1^2 + ...

    # Obtain the log-rate-density from the model:
    log_dN = model_log_dN(my_parameter_1, derived_quantity, ...)

    logl, log_norm, genq = model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    Turing.@addlogprob! logl
    Turing.@addlogprob! log_norm

    # Compute anything eles you need, maybe some more generated quantities
    additional_genq = (derived_quatity = derived_quantity, ...)

    # Return generated quantities 
    return merge(genq, additional_genq) # or just return genq if you have no additional generated quantities
"""
function model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
  Nobs = length(m1s)

  log_wts = map(m1s, m2s, log_wts) do mm1, mm2, log_wwt
      map(mm1, mm2, log_wwt) do m1, m2, log_wt
          log_dN(m1, m2) - log_wt
      end
  end

  Neff_samps = map(log_wts) do lws
      exp(2*logsumexp(lws) - logsumexp(2 .* lws))
  end

  log_likes = map(log_wts) do lw
      logsumexp(lw) - log(length(lw))
  end

  log_like_sum = sum(log_likes)

  log_wts_sel = map(m1s_sel, m2s_sel, log_pdraw) do m1, m2, logp
      log_dN(m1, m2) - logp
  end
  log_mu = logsumexp(log_wts_sel) - log(Ndraw)
  log_s2 = logsubexp(logsumexp(2 .* log_wts_sel) - 2.0*log(Ndraw), 2*log_mu - log(Ndraw))
  Neff_sel_est = exp(2*log_mu - log_s2)
  Neff_sel = 1/(1/Neff_sel_est + 1/Ndraw)

  log_norm_sum = -Nobs*log_mu

  mu = exp(log_mu)
  R = rand(Normal(Nobs/mu, sqrt(Nobs)/mu))
  
  m1_m2_draw = map(m1s, m2s, log_wts) do mm1, mm2, log_wwt
      i = rand(Categorical(exp.(log_wwt .- logsumexp(log_wwt))))
      (mm1[i], mm2[i])
  end

  m1s_draw = map(m1_m2_draw) do m1_m2
      m1_m2[1]
  end
  m2s_draw = map(m1_m2_draw) do m1_m2
      m1_m2[2]
  end

  i_sel = rand(Categorical(exp.(log_wts_sel .- logsumexp(log_wts_sel))))
  m1_draw_sel = m1s_sel[i_sel]
  m2_draw_sel = m2s_sel[i_sel]
  
  return log_like_sum, log_norm_sum, (R = R, Neff_sel = Neff_sel, Neff_samps = Neff_samps, m1s_popwt = m1s_draw, m2s_popwt = m2s_draw, m1_draw=m1_draw_sel, m2_draw=m2_draw_sel)
end

"""
  broken_pl_plp_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)

Returns a Turing model for our broken power law mass function, with power law
pairing function.

# Arguments

- `m1s`: An iterable containing iterables of the posterior samples for the
source-frame primary mass for our observations.
- `m2s`: Same, but for secondary mass.
- `log_wts`: Iterable with iterables of the log of the prior weight assigned to
each sample.
- `m1s_sel`: Samples of primary mass from the detected injections used to
estimate detector sensitivity.
- `m2s_sel`: Secondary mass.
- `log_pdraw`: The log of the (normalized) density from which the injections
were drawn.
- `Ndraw`: The number of injections drawn.

# Parameters

- `a1`: The low-mass power law slope.
- `a2`: The high-mass power law slope.
- `mb`: The break mass.
- `beta`: The pairing power law exponent.
- `R`: (generated) The volumetric merger rate per natural log mass squared at
`m1 = m2 = mb`.
- `Neff_sel`: (generated) The effective number of samples (see [Farr
(2019)](https://iopscience.iop.org/article/10.3847/2515-5172/ab1d5f)) in the
Monte-Carlo selection integral.
- `Neff_samps`: (generated) The effective number of samples in the posterior
samples after re-weighting to the current population model.
- `m1s_popwt`: (generated) The primary mass samples drawn from the population
model.
- `m2s_popwt`: (generated) The secondary mass samples drawn from the population
model.
- `m1_draw`: (generated) A draw of primary mass from the detected population.
- `m2_draw`: (generated) A draw of secondary mass from the detected population.
"""
@model function broken_pl_plp_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
  a1 ~ Uniform(0, 15)
  a2 ~ Uniform(-15, 0)
  mb ~ Uniform(sqrt(30), sqrt(200)) # Peaks somewhere around 10.

  beta ~ Uniform(-5, 5)

  log_r_ns ~ Uniform(log(0.01), log(100))
  r_ns = exp(log_r_ns)
  mu_ns ~ Uniform(1, 3)
  sigma_ns ~ Uniform(0.25, 2)

  log_dN = make_log_dN(BrokenPowerLaw(), PowerLawPairing(), a1, a2, mb, beta, r_ns, mu_ns, sigma_ns)

  logl_sum, lognorm_sum, generated_quantities = model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
  Turing.@addlogprob! logl_sum
  Turing.@addlogprob! lognorm_sum
  merge(generated_quantities, (r_ns=r_ns,))
end

@model function power_law_plus_gaussian_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
  a1 ~ Uniform(-10, 10)
  a2 ~ Uniform(-10, 10)
  mu ~ Uniform(5, 18)
  sigma ~ Uniform(0.25, 10)
  fg ~ Uniform(0, 1)
  beta ~ Uniform(-5, 5)

  log_r_ns ~ Uniform(log(0.01), log(100))
  r_ns = exp(log_r_ns)
  mu_ns ~ Uniform(1, 3)
  sigma_ns ~ Uniform(0.25, 2)

  log_dN = make_log_dN(PowerLawGaussian(), PowerLawPairing(), a1, a2, mu, sigma, fg, beta, r_ns, mu_ns, sigma_ns)

  logl_sum, lognorm_sum, generated_quantities = model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
  Turing.@addlogprob! logl_sum
  Turing.@addlogprob! lognorm_sum
  merge(generated_quantities, (r_ns=r_ns,))
end