const c = cosmology()

"""
    li_prior_wt(m1, m2, z)

Returns the LALInference prior in m1, m2, and z.

```math
p\\left( m_1, m_2, z \\right) \\propto \\left( 1 + z \\right)^2 d_L^2(z) \\frac{\\partial d_L}{\\partial z}
```
"""
function li_prior_wt(m1, m2, z)
    dL = ustrip(u"Gpc", luminosity_dist(c, z))
    dC = dL / (1+z)
    dH = ustrip(u"Gpc", hubble_dist(c, 0))
    ddLdz = (dC + (1+z)*dH/Cosmology.E(c, z))

    (1 + z)*(1 + z)*dL*dL*ddLdz
end

"""
    md_sfr(z[, lambda, z_p, kappa])

Returns the Madau-Dickinson SFR.

```math
\\frac{\\mathrm{d} N}{\\mathrm{d} V \\mathrm{d} t} \\propto \\frac{\\left( 1 + z \\right)^\\lambda}{1 + \\left( \\frac{1+z}{1+z_p}\\right)^\\kappa}
```

The default values are `(lambda, z_p, kappa) = (2.7, 1.9, 5.6)`.
"""
function md_sfr(z, lambda=2.7, z_p=1.9, kappa=5.6)
    zp1 = 1 + z
    zp1^lambda / (1 + (zp1 / (1 + z_p))^kappa)
end

"""
    md_sfr_zwt(z[, lambda, z_p, kappa])

Returns the weight in `z` corresponding to the Madau-Dickinson SFR.

```math
\\frac{\\mathrm{d} N}{\\mathrm{d} V \\mathrm{d} t} \\propto \\frac{\\left( 1 + z \\right)^\\lambda}{1 + \\left( \\frac{1+z}{1+z_p}\\right)^\\kappa}
```

The default values are `(lambda, z_p, kappa) = (2.7, 1.9, 5.6)`.
"""
function md_sfr_zwt(z, lambda=2.7, z_p=1.9, kappa=5.6)
    zp1 = 1 + z
    wt = zp1^lambda / (1 + (zp1 / (1 + z_p))^kappa)

    4*pi*md_sfr(z, lambda, z_p, kappa) * ustrip(u"Gpc^3", comoving_volume_element(c, z)) / zp1
end

"""
    read_selection(file)

Read the injection set from `file` and return `(m1s, m2s, zs, pdraw, Ndraw)` for sensitivity estimation.
"""
function read_selection(file)
    h5open(file, "r") do f
        i = f["injections"]
        m1s = convert(Vector{Float64}, read(i, "mass1_source"))
        m2s = convert(Vector{Float64}, read(i, "mass2_source"))
        zs = convert(Vector{Float64}, read(i, "redshift"))
    
        T = read(attributes(f), "analysis_time_s") / (356.25*24.0*3600.0) # yr

        pm1m2z = convert(Vector{Float64}, read(i, "mass1_source_mass2_source_sampling_pdf") .* read(i, "redshift_sampling_pdf"))
        pdraw = pm1m2z / T

        sel_flag = (read(i, "ifar_cwb") .> 1) .| (read(i, "ifar_gstlal") .> 1) .| (read(i, "ifar_mbta") .> 1) .| (read(i, "ifar_pycbc_bbh") .> 1) .| (read(i, "ifar_pycbc_hyperbank") .> 1)
    
        Ndraw = round(Int, read(attributes(i)["total_generated"]))
    
        m1s[sel_flag], m2s[sel_flag], zs[sel_flag], pdraw[sel_flag], Ndraw
    end
end

"""
    resample_selection(log_dN, m1s, m2s, pdraws, Ndraw)

Returns re-sampled `(m1s, m2s, pdraws, Ndraw)` as if injections had been drawn
from the given fiducial population.
"""
function resample_selection(log_dN, m1s, m2s, pdraws, Ndraw)
    pop_wts = exp.(log_dN.(m1s, m2s))
    wts = pop_wts ./ pdraws
    sum_wts = sum(wts)
    sum_wts2 = sum(wts.*wts)
    norm = sum_wts / Ndraw
    Neff_norm = round(Int, sum_wts*sum_wts / sum_wts2)

    inds = sample(1:length(m1s), weights(wts), Neff_norm)
    norm_pop_wts = pop_wts ./ norm
    m1s[inds], m2s[inds], norm_pop_wts[inds], Neff_norm
end

