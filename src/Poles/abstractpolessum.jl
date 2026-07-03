"""
    AbstractPolesSum

Supertype which represents a (block) function on the real axis as a sum of poles.
"""
abstract type AbstractPolesSum <: AbstractPoles end

amplitudes(P::AbstractPolesSum) = map(i -> amplitude(P, i), eachindex(P))

"""
    filling(P::AbstractPolesSum, μ::Real=0)

Calculate

```math
-\\frac{1}{π} \\int_{-∞}^μ \\mathrm{d}ω~\\mathrm{Im} P(ω+\\mathrm{i}0^+) .
```

If a pole exists at exactly ``μ``, its weight is counted half.
"""
function filling end

"""
    flip_spectrum!(P::AbstractPolesSum)

Reverse `P` and flip the sign of `locations(P)`.

See also [`flip_spectrum`](@ref).
"""
function flip_spectrum!(P::AbstractPolesSum)
    reverse!(P)
    l = locations(P)
    @. l = -l
    return P
end

"""
    flip_spectrum(P::AbstractPolesSum)

Reverse `P` and flip the sign of `locations(P)`.

See also [`flip_spectrum!`](@ref).
"""
flip_spectrum(P::AbstractPolesSum) = flip_spectrum!(copy(P))

"""
    merge_degenerate_poles!(P::AbstractPolesSum, tol::Real=0)

Merge poles whose locations are `≤ tol` apart.
"""
function merge_degenerate_poles! end

"""
    merge_negative_locations_to_zero!(P::AbstractPolesSum)

Find all `locations(P) <= 0` and merge them.
"""
function merge_negative_locations_to_zero!(P::AbstractPolesSum)
    # get information from P
    loc = locations(P)
    wgt = weights(P)
    idx_zeros = findall(<=(0), loc)
    isempty(idx_zeros) && return P
    # add up all weights
    w0 = sum(wgt[idx_zeros])
    i0 = popfirst!(idx_zeros)
    loc[i0] = 0
    wgt[i0] = w0
    # delete degenerate locations
    for i in reverse!(idx_zeros)
        deleteat!(loc, i)
        deleteat!(wgt, i)
    end
    return P
end

"""
    merge_small_weight!(P::AbstractPolesSum, tol::Real)

Merge poles with weight `<= tol` to its neighbors.

A given pole is split locally using the law of levers.
This conserves the zeroth and first moment for scalars.
"""
function merge_small_weight! end

"""
    moment(P::AbstractPolesSum, n::Int=0)

Return the `n`-th moment.
"""
function moment end

"""
    moments(P::AbstractPolesSum, ns)

Return the `n`-th moment for each `n` in `ns`.
"""
function moments(P::AbstractPolesSum, ns)
    return map(i -> moment(P, i), ns)
end

"""
    remove_zero_weight!(P::AbstractPolesSum, remove_zero::Bool=true)

Remove all poles which have zero weight.

If `remove_zero`, the pole at ``a_i = 0`` with zero weight is also removed.

See also [`remove_zero_weight`](@ref).
"""
function remove_zero_weight!(P::AbstractPolesSum, remove_zero::Bool = true)
    i = 1
    while i <= length(P)
        if iszero(location(P, i)) && !remove_zero
            # keep pole at origin
            i += 1
            continue
        end

        if iszero(weight(P, i))::Bool
            deleteat!(locations(P), i)
            deleteat!(weights(P), i)
        else
            i += 1
        end
    end
    return P
end

"""
    remove_zero_weight(P::AbstractPolesSum, remove_zero::Bool=true)

Remove all poles which have zero weight.

If `remove_zero`, the pole at ``a_i = 0`` with zero weight is also removed.

See also [`remove_zero_weight!`](@ref).
"""
function remove_zero_weight(P::AbstractPolesSum, remove_zero::Bool = true)
    return remove_zero_weight!(copy(P), remove_zero)
end

"""
    to_grid(P::AbstractPolesSum, grid::AbstractVector{<:Real})

Create a new [`AbstractPolesSum`](@ref) from `P` with locations given by `grid`.

A given pole is split locally conserving the zeroth and first moment.
If a pole is outside of `grid`, only the zeroth moment is conserved.
"""
function to_grid(P::AbstractPolesSum, grid::AbstractVector{<:Real})
    # check input
    _issorted_and_unique(grid)

    # new location and weights
    weights_new = [zero(first(weights(P))) for _ in eachindex(grid)]

    # run through each existing pole and split weight to new locations
    @inbounds for i in eachindex(P)
        loc = location(P, i)
        w = weight(P, i)
        if loc <= first(grid)
            # no pole to the left
            weights_new[begin] += w
        elseif loc >= last(grid)
            # no pole to the right
            weights_new[end] += w
        else
            # find next pole with higher location
            i = searchsortedfirst(grid, loc)
            if loc - grid[i - 1] < 10 * eps()
                # previous pole has same location
                weights_new[i - 1] += w
            elseif grid[i] - loc < 10 * eps()
                # current pole has same location
                weights_new[i] += w
            else
                # split such that zeroth and first moment is conserved
                loc_low = grid[i - 1]
                loc_high = grid[i]
                weights_new[i - 1] += (loc_high - loc) / (loc_high - loc_low) * w
                weights_new[i] += (loc - loc_low) / (loc_high - loc_low) * w
            end
        end
    end
    return typeof(P)(copy(grid), weights_new)
end

weight(P::AbstractPolesSum, i::Integer) = weights(P)[i]

weights(P::AbstractPolesSum) = P.weights

Base.eachindex(P::AbstractPolesSum) = eachindex(locations(P))

function Base.iterate(P::AbstractPolesSum, i = 0)
    @inline
    next = i + 1
    (i == length(P)) && return nothing
    return ((location(P, next), weight(P, next)), next)
end

Base.reverse(P::AbstractPolesSum) = reverse!(copy(P))

function Base.reverse!(P::AbstractPolesSum)
    reverse!(locations(P))
    reverse!(weights(P))
    return P
end

function Base.sort!(P::AbstractPolesSum)
    p = sortperm(locations(P))
    P.locations[:] = P.locations[p]
    P.weights[:] = P.weights[p]
    return P
end

function LinearAlgebra.rmul!(P::AbstractPolesSum, α::Number)
    rmul!(weights(P), α::Number)
    return P
end
