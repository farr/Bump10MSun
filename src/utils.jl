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

function cumtrapz(xs, ys)
    dx = xs[2:end] - xs[1:end-1]
    integral = zeros(typeof(xs[1]), length(xs))
    integral[2:end] = cumsum((ys[1:end-1] .+ ys[2:end]) .* dx) / 2
    integral
end
