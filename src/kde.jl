# Required because KernelDensity.jl uses axis-aligned KDEs.


struct KDE <: ContinuousMultivariateDistribution
    pts::Matrix{Float64}
    chol_bw::Cholesky{Float64, Matrix{Float64}}
end

"""
    KDE(pts)

Return a kernel density object out of the given points.

`pts` should have shape `(ndim, npts)` or `(npts,)` for a 1-dimensional KDE.
"""
function KDE(pts::Matrix{Float64})
    nd, np = size(pts)
    S = cov(pts') ./ np^(2/(nd+4))

    KDE(pts, cholesky(S))
end

function KDE(pts::Vector{Float64})
    KDE(reshape(pts, 1, :))
end

"""The number of dimensions in the KDE."""
function ndim(k::KDE)
    size(k.pts, 1)
end
"""The number of points stored in the KDE."""
function npts(k::KDE)
    size(k.pts, 2)
end
"""The cholesky factor of the KDE bandwidth matrix."""
function chol_bw(k::KDE)
    k.chol_bw
end

Distributions.length(k::KDE) = ndim(k)
Distributions.sampler(k::KDE) = k
Distributions.eltype(::KDE) = Float64
function Distributions._rand!(rng::AbstractRNG, k::KDE, x::AbstractVector)
    i = rand(rng, 1:npts(k))
    x .= k.pts[:, i] .+ k.chol_bw.L * randn(rng, ndim(k))
end
function Distributions._logpdf(k::KDE, x::AbstractArray)
    lp = -Inf
    for j in 1:npts(k)
        r = x .- k.pts[:, j]
        lp = logaddexp(lp, -(r' * (k.chol_bw \ r))/2)
    end

    lp = lp - sum(log.(sqrt(2*pi).*diag(k.chol_bw.L))) - log(npts(k))
    lp
end

"""
    logpdf_credible_levels(k, pts, qs)

Return the logpdf values of the KDE `k` corresponding to the credible levels
`qs` of the density of points `pts`.
"""
function logpdf_credible_levels(k::KDE, pts::Matrix{Float64}, qs::Vector{Float64})
    np = size(pts, 2)

    lpdf = [Distributions.logpdf(k, pts[:, i]) for i in 1:np]
    quantile.((lpdf,), qs)
end