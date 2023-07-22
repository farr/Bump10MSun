"""    with_seed(f, seed)

Perform the zero-argument function `f` after seeding the default RNG with `seed`.

Ensures that the RNG is re-seeded randomly when `f` exits normally or abnormally
(i.e. by throwing an exception).

Note the useful syntax:

    with_seed(1234) do
        ...
    end
"""
function with_seed(f, seed)
    Random.seed!(seed)
    try
        f()
    finally
        Random.seed!()
    end
end

"""
    cumtrapz(xs, ys)

Return an array of the cumulative trapezoidal integral of `ys` with respect to
`xs`.

The first element of the array will always be `0`.
"""
function cumtrapz(xs, ys)
    dx = xs[2:end] - xs[1:end-1]
    integral = zeros(typeof(xs[1]), length(xs))
    integral[2:end] = cumsum((ys[1:end-1] .+ ys[2:end]) .* dx) / 2
    integral
end

"""
    load_pe_samples(dirname, glob_pattern, dataset_key)

Return `(samples, filenames, gwnames)` from all HDF5 files matching
`glob_pattern` in `dirname`.

The array `samples` contains the dataset from each HDF5 file corresponding to
dataset_key.  `gwnames` is an array of names of the form `GWYYMMDD[_NNNNNN]` for
each file.

# Examples
```julia
samples, fnames, gwnames = load_pe_samples("/path/to/PE/directory", "GW*.h5", "PublicationSamples/posterior_samples")
```
"""
function load_pe_samples(dirname, glob_pattern, dataset_key)
    samps = []
    fnames = []
    for fn in glob(glob_pattern, dirname)
        h5open(fn, "r") do f
            push!(samps, read(f, dataset_key))
            push!(fnames, fn)
        end
    end
    gwnames = [String(match(r"^.*(GW[0-9]+[_]*[0-9]*).*$", f)[1]) for f in fnames]

    return samps, fnames, gwnames
end

"""
    chirp_mass(m1, m2)

Return the chirp mass associated with comp
"""
function chirp_mass(m1, m2)
    mt = m1+m2
    eta = m1*m2/(mt*mt)

    mt*eta^(3/5)
end