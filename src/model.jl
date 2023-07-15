"""The minimum mass considered."""
const mlow = 3.0

"""The maximum mass considered."""
const mhigh = 20.0

"""The reference redshift, at which rates are evaluated."""
const zref = 0.0

function _fm(m, a1, a2, mb)
    x = m/mb
    if m < mlow || m > mhigh
        -Inf
    elseif m < mb
        a1*log(x)
    else
        a2*log(x)
    end
end

function _fm(m, a1, a2, a3, mb12, mb23)
    if m < mlow || m > mhigh
        -Inf
    elseif m < mb12
        a1*log(m/mb12)
    elseif m < mb23
        a2*log(m/mb12)
    else
        a3*log(m/mb12) + (a2-a3)*log(mb23/mb12)
    end
end

"""
    make_log_dN(a1, a2, mb, mu_q, sigma_q)

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
"""
square(x) = x*x
function make_log_dN(a1, a2, mb, mu_q, sigma_q)
    function log_dN(m1, m2)
        q = m2/m1
        _fm(m1, a1, a2, mb) + _fm(m2, a1, a2, mb) - log(m1) - log(m2) - 0.5*square((q-mu_q)/sigma_q) + 0.5*square((1-mu_q)/sigma_q)
    end
    log_dN
end

"""
    make_log_dN(a1, a2, mb, beta)

Produces a function that computes `(m1, m2) -> log_dNdm1dm2` for the
single-break power law population model.

The population model is 

```math
m_1 m_2 \\frac{\\mathrm{d}N}{\\mathrm{d} m_1 \\mathrm{d} m_2 \\mathrm{d} V \\mathrm{d} t} = R f_m\\left( m_1 \\right) f_m\\left( m_2 \\right) q^\\beta
```
where
```math
f_m\\left( m \\right) = \\begin{cases}
m^{\\alpha_1} & m_\\mathrm{min} < m < m_b \\
m^{\\alpha_2} & m_b \\leq m < m_\\mathrm{max}
\\end{cases}
```
is the "common" black hole mass function, and a power-law pairing function in
``q`` is used with parameter ``beta``.  The function returned is normalized so
that ``R`` measures the volumetric merger rate per natural log mass squared at
`m1 = m2 = m_b`.
"""
function make_log_dN(a1, a2, mb, beta)
    function log_dN(m1, m2)
        q = m2/m1
        _fm(m1, a1, a2, mb) + _fm(m2, a1, a2, mb) - log(m1) - log(m2) + beta*log(q)
    end
    log_dN
end

"""
    make_log_dN(a1, a2, a3, mb12, mb23, mu_q, sigma_q)

Produces a function that computes `(m1, m2) -> log_dNdm1dm2` for the two-break
power law population model.

The population model is 

```math
m_1 m_2 \\frac{\\mathrm{d}N}{\\mathrm{d} m_1 \\mathrm{d} m_2 \\mathrm{d} V \\mathrm{d} t} = R f_m\\left( m_1 \\right) f_m\\left( m_2 \\right) \\exp\\left( -\\frac{\\left(q - \\mu_q\\right)^2}{2 \\sigma_q^2} + \\frac{\\left(1-\\mu_q\\right)^2}{2 \\sigma_q^2} \\right)^2
```
where
```math
f_m\\left( m \\right) = \\begin{cases}
m^{\\alpha_1} & m_\\mathrm{min} < m < m_{b,12} \\
m^{\\alpha_2} & m_{b,12} < m < m_{b,23} \\
m^{\\alpha_3} & m_{b,23} \\leq m < m_\\mathrm{max}
\\end{cases}
```
is the "common" black hole mass function, and a Gaussian-shaped pairing function
in ``q`` is used with parameters ``mu_q`` and ``sigma_q``.  The function
returned is normalized so that ``R`` measures the volumetric merger rate per
natural log mass squared at `m1 = m2 = mb12`.
"""
function make_log_dN(a1, a2, a3, mb12, mb23, mu_q, sigma_q)
    function log_dN(m1, m2)
        q = m2/m1
        _fm(m1, a1, a2, a3, mb12, mb23) + _fm(m2, a1, a2, a3, mb12, mb23) - log(m1) - log(m2) - 0.5*square((q-mu_q)/sigma_q) + 0.5*square((1-mu_q)/sigma_q)
    end
    log_dN
end

"""
    make_log_dN(a1, a2, a3, mb12, mb23, beta)

Produces a function that computes `(m1, m2) -> log_dNdm1dm2` for the two-break
power law population model.

The population model is 

```math
m_1 m_2 \\frac{\\mathrm{d}N}{\\mathrm{d} m_1 \\mathrm{d} m_2 \\mathrm{d} V \\mathrm{d} t} = R f_m\\left( m_1 \\right) f_m\\left( m_2 \\right) q^\\beta
```
where
```math
f_m\\left( m \\right) = \\begin{cases}
m^{\\alpha_1} & m_\\mathrm{min} < m < m_{b,12} \\
m^{\\alpha_2} & m_{b,12} < m < m_{b,23} \\
m^{\\alpha_3} & m_{b,23} \\leq m < m_\\mathrm{max}
\\end{cases}
```
is the "common" black hole mass function, and a power-law pairing function in
``q`` with parameter ``beta``.  The function returned is normalized so that
``R`` measures the volumetric merger rate per natural log mass squared at `m1 =
m2 = mb12`.
"""
function make_log_dN(a1, a2, a3, mb12, mb23, beta)
    function log_dN(m1, m2)
        q = m2/m1
        _fm(m1, a1, a2, a3, mb12, mb23) + _fm(m2, a1, a2, a3, mb12, mb23) - log(m1) - log(m2) + beta*log(q)
    end
    log_dN
end

"""
    make_dNdm1(args...)

Returns a function that computes the marginal density (unnormalized) in `m1`.
Arguments are passed to ``make_log_dN``.  
"""
function make_dNdm1(arg1, args...)
    log_dN = make_log_dN(arg1, args...)
    make_dNdm1(log_dN)
end
function make_dNdm1(log_dN)
    m1 -> begin
        if m1 <= mlow
            zero(m1)
        elseif m1 > mhigh
            zero(m1)
        else
            n = Int(round(100*(log(m1)-log(mlow)))) + 2
            ms = exp.(range(log(mlow), stop=log(m1), length=n))
            dN = exp.(map(m2->log_dN(m1,m2), ms))
            trapz(ms, dN)
        end
    end
end

"""
    make_dNdm2(args...)

Returns a function that computes the marginal density (unnormalized) in `m2`;
``args`` are passed to ``make_log_dN``.
"""
function make_dNdm2(arg1, args...)
    log_dN = make_log_dN(arg1, args...)
    make_dNdm2(log_dN)
end
function make_dNdm2(log_dN)
    m2 -> begin
        if m2 < mlow
            zero(m2)
        elseif m2 > mhigh
            zero(m2)
        else
            n = Int(round(100*(log(mhigh) - log(m2)))) + 2
            ms = exp.(range(log(m2), stop=log(mhigh), length=n))
            dN = exp.(map(m1->log_dN(m1, m2), ms))
            trapz(ms, dN)
        end
    end
end

"""
Returns a function that computes the marginal (unnormalized) density in `q`.
"""
function make_dNdq(arg1, args...)
    log_dN = make_log_dN(arg1, args...)
    make_dNdq(log_dN)
end
function make_dNdq(log_dN)
    q -> begin
        if q < mlow/mhigh
            zero(q)
        elseif q > 1
            zero(q)
        else
            ml = mlow/q
            mh = mhigh
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
function make_dNdm(a1, a2, mb)
    m -> begin
        exp(_fm(m, a1, a2, mb))/m     
    end
end
function make_dNdm(a1, a2, a3, mb12, mb23)
    m -> begin
        exp(_fm(m, a1, a2, a3, mb12, mb23))/m
    end
end

"""
Returns the "pairing probability", the function ``g(q)``, either a Gaussian with
``mu_q`` and ``sigma_q`` parameters, or a power law with parameter ``beta``.
"""
function make_pairing_prob(mu_q, sigma_q)
    q -> begin
        exp(-0.5*square((q-mu_q)/sigma_q) + 0.5*square((1-mu_q)/sigma_q))
    end
end

function make_pairing_prob(beta)
    q -> begin
        q^beta
    end 
end

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
    Neff_sel = exp(2*log_mu - log_s2)

    log_norm_sum = -Nobs*log_mu

    mu_R = Nobs/exp(log_mu)
    sigma_R = sqrt(Nobs) / exp(log_mu)

    R = rand(Normal(mu_R, sigma_R))
    
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
    broken_pl_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)

Returns a Turing model for our broken power law mass function with a Gaussian
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
- `mu_q`: The peak of the pairing function Gaussian.
- `sigma_q`: The width of the pairing function Gaussian.
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
@model function broken_pl_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    a1 ~ Uniform(0, 15)
    a2 ~ Uniform(-15, 0)
    mb ~ Uniform(sqrt(30), sqrt(200)) # Peaks somewhere around 10.
    mu_q ~ Uniform(0.2, 1.0)
    sigma_q ~ Uniform(0.1, 2)

    log_dN = make_log_dN(a1, a2, mb, mu_q, sigma_q)

    logl_sum, lognorm_sum, generated_quantities = model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    Turing.@addlogprob! logl_sum
    Turing.@addlogprob! lognorm_sum
    generated_quantities
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

    log_dN = make_log_dN(a1, a2, mb, beta)

    logl_sum, lognorm_sum, generated_quantities = model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    Turing.@addlogprob! logl_sum
    Turing.@addlogprob! lognorm_sum
    generated_quantities
end

"""
    two_broken_pl_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)

Returns a Turing model for our broken power law mass function with a Gaussian
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
- `a2`: The mid-mass power law slope.
- `a3`: The high-mass power law slope.
- `mb12`: The low-mid break mass.
- `mb23`: The mid-high break mass.
- `mu_q`: The peak of the pairing function Gaussian.
- `sigma_q`: The width of the pairing function Gaussian.
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
@model function two_broken_pl_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    a1 ~ Uniform(-10, 10)
    a2 ~ Uniform(-10, 10)
    a3 ~ Uniform(-10, 10)

    # Flat in mb12-mb23 space.
    x12 ~ Uniform(0,1)
    x23 ~ Uniform(0,1)
    mb12 = mlow+1 + (mhigh-1 - (mlow+1))*x12
    mb23 = mb12 + (mhigh-1 - mb12)*x23
    Turing.@addlogprob! log((mhigh-1) - (mlow+1)) + log(mhigh-1 - mb12)

    mu_q ~ Uniform(0.2, 1.0)
    sigma_q ~ Uniform(0.1, 2)

    log_dN = make_log_dN(a1, a2, a3, mb12, mb23, mu_q, sigma_q)

    logl_sum, lognorm_sum, generated_quantities = model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    Turing.@addlogprob! logl_sum
    Turing.@addlogprob! lognorm_sum
    merge(generated_quantities, (mb12=mb12, mb23=mb23))
end

"""
    two_broken_pl_plp_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)

Returns a Turing model for our broken power law mass function with a power law
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
- `a2`: The mid-mass power law slope.
- `a3`: The high-mass power law slope.
- `mb12`: The low-mid break mass.
- `mb23`: The mid-high break mass.
- `beta`: The power law index of the pairing function.
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
@model function two_broken_pl_plp_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    a1 ~ Uniform(-10, 10)
    a2 ~ Uniform(-10, 10)
    a3 ~ Uniform(-10, 10)

    # Flat in mb12-mb23 space.
    x12 ~ Uniform(0,1)
    x23 ~ Uniform(0,1)
    mb12 = mlow+1 + (mhigh-1 - (mlow+1))*x12
    mb23 = mb12 + (mhigh-1 - mb12)*x23
    Turing.@addlogprob! log((mhigh-1) - (mlow+1)) + log(mhigh-1 - mb12)

    beta ~ Uniform(-5, 5)

    log_dN = make_log_dN(a1, a2, a3, mb12, mb23, beta)

    logl_sum, lognorm_sum, generated_quantities = model_body(log_dN, m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    Turing.@addlogprob! logl_sum
    Turing.@addlogprob! lognorm_sum
    merge(generated_quantities, (mb12=mb12, mb23=mb23))
end
