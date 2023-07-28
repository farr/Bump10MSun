"""
    logpdf_credible_levels(k, pts, qs)

Return the logpdf values of the KDE `k` corresponding to the credible levels
`qs` of the density of points `pts`.
"""
function logpdf_credible_levels(k, pts, qs)
    np = size(pts, 2)

    lpdf = [Distributions.logpdf(k, pts[:, i]) for i in 1:np]
    quantile.((lpdf,), qs)
end