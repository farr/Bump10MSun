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
    
        s1x = convert(Vector{Float64}, read(i, "spin1x"))
        s1y = convert(Vector{Float64}, read(i, "spin1y"))
        s1z = convert(Vector{Float64}, read(i, "spin1z"))

        s2x = convert(Vector{Float64}, read(i, "spin2x"))
        s2y = convert(Vector{Float64}, read(i, "spin2y"))
        s2z = convert(Vector{Float64}, read(i, "spin2z"))

        a1 = sqrt.(s1x .* s1x .+ s1y .* s1y .+ s1z .* s1z)
        a2 = sqrt.(s2x .* s2x .+ s2y .* s2y .+ s2z .* s2z)

        # p(a, Omega) = 1/(4*pi) 
        # p(xyz) d(xyz) = p(a, Omega) da dOmega
        # p(xyz) a^2 da dOmega = p(a, Omega) da dOmega
        # p(xyz) = p(a, Omega) / a^2

        pxyz1 = 1 ./ (4 .* pi .* a1 .* a1)
        pxyz2 = 1 ./ (4 .* pi .* a2 .* a2)

        T = read(attributes(f), "analysis_time_s") / (356.25*24.0*3600.0) # yr

        pm1m2zs1s2 = convert(Vector{Float64}, read(i, "sampling_pdf"))
        # p(m1, m2, z, s1, s1) / p(s1) / p(s2) / T => p(m1, m2, z, t)
        pdraw = pm1m2zs1s2 ./ T ./ pxyz1 ./ pxyz2

        sel_flag_o3 = (read(i, "name") == "o3") .& ((read(i, "ifar_cwb") .> 1/far_threshold) .| (read(i, "ifar_gstlal") .> 1/far_threshold) .| (read(i, "ifar_mbta") .> 1/far_threshold) .| (read(i, "ifar_pycbc_bbh") .> 1/far_threshold) .| (read(i, "ifar_pycbc_hyperbank") .> 1/far_threshold))
        sel_flag_o1o2 = (read(i, "name") !== "o3") .& (read(i, "optimal_snr_net") .> snr_threshold)
        sel_flag = sel_flag_o3 .| sel_flag_o1o2
    
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

