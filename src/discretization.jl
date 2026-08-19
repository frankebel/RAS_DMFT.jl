"""
    discretize_to_grid(
        f::AbstractVector{<:Real}, W::AbstractVector{<:Real}, grid::AbstractVector{<:R}
    ) where {R<:Real}

Discretize the given function `f` defined on locations `W` to `grid`.

Returns a [`PolesSum`](@ref) object with `locations(P) == grid`.
"""
function discretize_to_grid(
        f::AbstractVector{<:Real}, W::AbstractVector{<:Real}, grid::AbstractVector{<:R}
    ) where {R <: Real}
    # check input
    _issorted_and_unique(grid)
    eachindex(f) == eachindex(W) || throw(ArgumentError("f and W must have same indexing"))
    Base.require_one_based_indexing(grid)

    weights = zero(grid)
    ig = firstindex(grid)
    border_next = (ig + 1) > lastindex(grid) ? typemax(R) : (grid[1] + grid[2]) / 2
    for i in eachindex(f)
        i == firstindex(f) && continue # no point to the left
        if W[i] > border_next
            # next location in grid
            ig += 1
            border_next = ig == lastindex(grid) ? typemax(R) : (grid[ig] + grid[ig + 1]) / 2
        end
        # add area by trapezoidal rule
        weights[ig] += (W[i] - W[i - 1]) * (f[i - 1] + f[i]) / 2
    end
    weights ./= π # correct norm

    return PolesSum(copy(grid), weights)
end
