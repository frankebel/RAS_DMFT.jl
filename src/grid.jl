# creation, manipulation of grids

"""
    grid_log(ω_max::Real, Λ::Real, n::Int)

Create `n` sorted poles on a logarithmic grid with highest value `ω_max`.

The poles have constant ratio

```math
\\frac{ω_i}{ω_{i+1}} =  \\frac{1}{Λ}.
```
"""
function grid_log(ω_max::Real, Λ::Real, n::Int)
    ω_max >= 0 || throw(ArgumentError("negative frequency ω_max"))
    Λ > 1 || throw(ArgumentError("invalid discretization parameter Λ"))
    n >= 1 || throw(ArgumentError("invalid number of points n"))
    result = map(i -> Λ^(-i) * ω_max, 0.0:(n - 1))
    reverse!(result) # from small to big values
    return result
end

"""
    grid_interpolate(a::AbstractVector{<:R}, n::Int) where {R<:Real}

For each point ``a_i`` linearly interpolate `n` point in the interval

```math
\\left[\\frac{a_{i-1} + a_i}{2}, \\frac{a_i + a_{i+1}}{2}\\right].
```
"""
function grid_interpolate(a::AbstractVector{<:Real}, n::Int)
    # check input
    _issorted_and_unique(a)
    Base.require_one_based_indexing(a)
    n >= 1 || throw(ArgumentError("n must be positive"))

    result = similar(a, length(a) * n + 1)
    for i in eachindex(a)
        # calculate end points of interval
        if i == firstindex(a)
            ω_low = a[i]
            ω_high = 0.5 * (a[i] + a[i + 1])
        elseif i == lastindex(a)
            ω_low = 0.5 * (a[i - 1] + a[i])
            ω_high = a[i]
        else
            ω_low = 0.5 * (a[i - 1] + a[i])
            ω_high = 0.5 * (a[i] + a[i + 1])
        end
        Δ = ω_high - ω_low # interval width
        # calculate points
        result[(1 + (i - 1) * n):(i * n)] = range(ω_low; step = Δ / n, length = n)
    end
    result[end] = last(a)

    return result
end
