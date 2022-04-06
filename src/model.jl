"""The minimum mass considered."""
const mlow = 1.0

"""The maximum mass considered."""
const mhigh = 20.0

"""The reference redshift, at which rates are evaluated."""
const zref = 0.4

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

"""
    make_log_dN(a1, a2, mb, beta)

Produces a function that computes `(m1, m2) -> log_dNdm1dm2` for the population
model.

The population model is 

```math
m_1 m_2 \\frac{\\mathrm{d}N}{\\mathrm{d} m_1 \\mathrm{d} m_2 \\mathrm{d} V \\mathrm{d} t} = R f_m\\left( m_1 \\right) f_m\\left( m_2 \\right) \\left(\\frac{m_2}{m_1}\\right)^\\beta
```
where
```math
f_m\\left( m \\right) = \\begin{cases}
m^{\\alpha_1} & m_\\mathrm{min} < m < m_b \\
m^{\\alpha_2} & m_b \\leq m < m_\\mathrm{max}
\\end{cases}
```
is the "common" black hole mass function, and the term involving ``\\beta`` is a
"pairing function" modification.  The function returned is normalized so that
``R`` measures the volumetric merger rate per natural log mass squared at `m1 =
m2 = m_b`.
"""
function make_log_dN(a1, a2, mb, beta)
    function log_dN(m1, m2)
        _fm(m1, a1, a2, mb) + _fm(m2, a1, a2, mb) - log(m1) - log(m2) + beta*log(m2/m1)
    end
    log_dN
end

function _I1(m1, a1, a2, mb, beta)
    if m1 < mlow
        zero(m1)
    elseif m1 < mb
        m1^(-beta)*mb^(-a1)*(m1^(a1+beta) - mlow^(a1+beta))/(a1+beta)
    else
        m1^(-beta)*mb^(-a1)*(mb^(a1+beta) - mlow^(a1+beta))/(a1+beta)
    end
end

function _I2(m1, a1, a2, mb, beta)
    if m1 < mb
        zero(m1)
    elseif m1 < mhigh
        ((m1/mb)^a2 - (mb/m1)^beta)/(a2+beta)
    else
        ((mhigh/mb)^a2 - (mb/mhigh)^beta)/(a2+beta)
    end
end

"""
    make_dNdm1(a1, a2, mb, beta)

Returns a function that computes the marginal density (unnormalized) in `m1`.
"""
function make_dNdm1(a1, a2, mb, beta)
    m1 -> begin
        I = _I1(m1, a1, a2, mb, beta) + _I2(m1, a1, a2, mb, beta)
        if m1 < mlow
            zero(m1)
        elseif m1 < mb
            I*(m1/mb)^a1/m1
        elseif m1 < mhigh
            I*(m1/mb)^a2/m1
        else
            zero(m1)
        end
    end
end

function _J1(m2, a1, a2, mb, beta)
    if m2 < mlow
        m2^beta*mb^(-a1)*(mb^(a1-beta) - mlow^(a1-beta))/(a1-beta)
    elseif m2 < mb
        ((m2/mb)^beta - (m2/mb)^a1)/(a1-beta)
    else
        zero(m2)
    end
end

function _J2(m2, a1, a2, mb, beta)
    if m2 < mb
        m2^beta*mb^(-a2)*(mhigh^(a2-beta)-mb^(a2-beta))/(a2-beta)
    elseif m2 < mhigh
        m2^beta*mb^(-a2)*(mhigh^(a2-beta) - m2^(a2-beta))/(a2-beta)
    else
        zero(m2)
    end
end

"""
    make_dNdm2(a1, a2, mb, beta)

Returns a function that computes the marginal density (unnormalized) in `m2`.
"""
function make_dNdm2(a1, a2, mb, beta)
    m2 -> begin
        I = _J1(m2, a1, a2, mb, beta) + _J2(m2, a1, a2, mb, beta)
        if m2 < mlow
            zero(m2)
        elseif m2 < mb
            I*(m2/mb)^a1/m2
        elseif m2 < mhigh
            I*(m2/mb)^a2/m2
        else
            zero(m2)
        end
    end
end

"""
    broken_pl_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)

Returns a Turing model for our broken power law mass function.

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
- `Neff_sel`: The effective number of samples (see [Farr
  (2019)](https://iopscience.iop.org/article/10.3847/2515-5172/ab1d5f)) in the
  Monte-Carlo selection integral.
"""
@model function broken_pl_model(m1s, m2s, log_wts, m1s_sel, m2s_sel, log_pdraw, Ndraw)
    Nobs = length(m1s)
    a1 ~ Uniform(0, 10)
    a2 ~ Uniform(-10, 0)
    mb ~ Uniform(sqrt(30), sqrt(200)) # Peaks somewhere around 10.
    beta ~ Uniform(-2, 6)

    log_dN = make_log_dN(a1, a2, mb, beta)

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

    Turing.@addlogprob! sum(log_likes)

    log_wts_sel = map(m1s_sel, m2s_sel, log_pdraw) do m1, m2, logp
        log_dN(m1, m2) - logp
    end
    log_mu = logsumexp(log_wts_sel) - log(Ndraw)
    log_s2 = logsubexp(logsumexp(2 .* log_wts_sel) - 2.0*log(Ndraw), 2*log_mu - log(Ndraw))
    Neff_sel = exp(2*log_mu - log_s2)

    Turing.@addlogprob! -Nobs*log_mu

    mu_R = Nobs/exp(log_mu)
    sigma_R = sqrt(Nobs) / exp(log_mu)

    R = rand(Normal(mu_R, sigma_R))
    
    return (R = R, Neff_sel = Neff_sel, Neff_samps = Neff_samps)
end