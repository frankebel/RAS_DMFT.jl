"""
    AbstractPolesSum{A, B} <: AbstractPoles{A, B}

Supertype which represents a (block) function on the real axis as a sum of poles.

The canonical representation requires the locations to be strictly increasing
(sorted with no duplicates).
"""
abstract type AbstractPolesSum{A, B} <: AbstractPoles{A, B} end

amplitudes(P::AbstractPolesSum, args...; kwargs...) = map(i -> amplitude(P, i, args...; kwargs...), eachindex(P))

"""
    anderson_matrix(P::AbstractPolesSum)

Calculate scaling factor ``B_0`` and Anderson matrix ``H_\\mathrm{A}``
given a sum of poles `P`.

Reference: [DOI](https://doi.org/10.48550/arXiv.2605.04974), appendix A3d
"""
function anderson_matrix end

"""
    arrowhead_matrix(P::AbstractPolesSum, args...; kwargs...)

Calculate the (block) arrowhead matrix representation of

```math
\\frac{1}{z - \\mathbb{0} - P(z)} .
```

See also [`amplitude`](@ref) for details of `args...` and `kwargs...`.

```jldoctest
julia> P = PolesSum(1:2, [4, 9])
PolesSum{Int64, Int64} with 2 poles

julia> arrowhead_matrix(P)
3×3 Matrix{Float64}:
 0.0  2.0  3.0
 2.0  1.0  0.0
 3.0  0.0  2.0
```
"""
function arrowhead_matrix end

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
function merge_degenerate_poles!(P::AbstractPolesSum, tol::Real = 0)
    # check input
    tol >= 0 || throw(ArgumentError("tol must not be negative"))
    # get information from P
    locs = locations(P)
    wgts = weights(P)
    # pole(s) at [-tol, tol]
    idx_zeros = findall(i -> abs(i) <= tol, locs)
    if !isempty(idx_zeros)
        i0 = popfirst!(idx_zeros)
        locs[i0] = 0
        for i in reverse!(idx_zeros)
            wgts[i0] = _axpy!(true, popat!(wgts, i), wgts[i0])
            deleteat!(locs, i)
        end
    end
    # pole(s) at tol → ∞
    i = findfirst(>(0), locs)
    isnothing(i) && (i = lastindex(locs)) # enforce `i` to be a number
    while i < lastindex(locs)
        if locs[i + 1] - locs[i] <= tol
            # merge
            wgts[i] = _axpy!(true, popat!(wgts, i + 1), wgts[i])
            deleteat!(locs, i + 1) # keep location closer to zero
        else
            # increment index
            i += 1
        end
    end
    # pole(s) at -tol → -∞
    i = findlast(<(0), locs)
    isnothing(i) && (i = firstindex(locs)) # enforce `i` to be a number
    while i > firstindex(locs)
        if locs[i] - locs[i - 1] <= tol
            # merge
            wgts[i - 1] = _axpy!(true, popat!(wgts, i), wgts[i - 1])
            deleteat!(locs, i - 1) # keep location closer to zero
            i -= 1
        else
            # decrement index
            i -= 1
        end
    end
    return P
end

# y ← α*x + y for both scalar and block weights.
_axpy!(α, x::Number, y::Number) = α * x + y
_axpy!(α, x::AbstractArray, y::AbstractArray) = axpy!(α, x, y)

"""
    merge_negative_locations_to_zero!(P::AbstractPolesSum)

Find all `locations(P) <= 0` and merge them.
"""
function merge_negative_locations_to_zero!(P::AbstractPolesSum)
    # get information from P
    locs = locations(P)
    wgts = weights(P)
    idx_zeros = findall(<=(0), locs)
    isempty(idx_zeros) && return P
    # add up all weights
    w0 = sum(wgts[idx_zeros])
    i0 = popfirst!(idx_zeros)
    locs[i0] = 0
    wgts[i0] = w0
    # delete degenerate locations
    for i in reverse!(idx_zeros)
        deleteat!(locs, i)
        deleteat!(wgts, i)
    end
    return P
end

"""
    merge_small_weight!(P::AbstractPolesSum, tol::Real)

Merge poles with weight `<= tol` to its neighbors.

A given pole is split locally using the law of levers.
This conserves the zeroth and first moment for scalars.
For block weights, the pole size is measured by the largest eigenvalue.
"""
function merge_small_weight!(P::AbstractPolesSum, tol::Real)
    # check input
    tol >= 0 || throw(ArgumentError("negative tol is invalid"))
    # loop over all poles
    i = 1
    while i <= length(P)
        loc = location(P, i)
        wgt = weight(P, i)
        if _mag(wgt) > tol
            # enough weight, go to next
            i += 1
            continue
        end
        if i == 1
            # add weight to next pole
            weights(P)[2] = _axpy!(true, wgt, weights(P)[2])
            deleteat!(locations(P), 1)
            deleteat!(weights(P), 1)
        elseif i == length(P)
            # add weight to previous pole
            weights(P)[end - 1] = _axpy!(true, wgt, weights(P)[end - 1])
            pop!(locations(P))
            pop!(weights(P))
        else
            # split weight such that zeroth and first moment is conserved
            loc_prev = location(P, i - 1)
            loc_next = location(P, i + 1)
            α = (loc_next - loc) / (loc_next - loc_prev)
            weights(P)[i - 1] = _axpy!(α, wgt, weights(P)[i - 1])
            weights(P)[i + 1] = _axpy!(1 - α, wgt, weights(P)[i + 1])
            deleteat!(locations(P), i)
            deleteat!(weights(P), i)
        end
    end
    return P
end

# Size of a (block) weight for the small-weight threshold.
_mag(x::Number) = x
_mag(x::AbstractArray) = eigmax(Hermitian(x))

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
    shift_spectrum!(P::AbstractPolesSum, μ::Real)

Shift locations of `P` by ``a_i → a_i - μ``.
"""
function shift_spectrum!(P::AbstractPolesSum, μ::Real)
    locations(P) .-= μ
    return P
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
    @inbounds for (loc, w) in P
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
