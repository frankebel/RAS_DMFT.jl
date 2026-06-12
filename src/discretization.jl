# discretize given object to `n` poles

"""
    discretize_similar_weight(P::PolesSum, δ0::Real, n::Int)

Discretize `P` to `n` poles such that each new pole has approximately equal weight.

For the new pole at location zero, all weights of `P` in ``[-δ_0, δ_0]`` are summed up.
"""
function discretize_similar_weight(P::PolesSum, δ0::Real, n::Integer)
    l_old = locations(P)
    w_old = weights(P)
    δ0 >= 0 || throw(ArgumentError("negative δ0"))
    n >= 3 || throw(ArgumentError("at least 3 poles necessary"))
    isodd(n) || throw(ArgumentError("number of poles must be odd"))

    T = eltype(P)
    loc_new = Vector{T}(undef, n)
    wgt_new = Vector{T}(undef, n)

    # weight at zero is weight in [-δ0, -δ0]
    i_minus = searchsortedfirst(l_old, -δ0)
    i_plus = searchsortedlast(l_old, δ0)
    w0 = sum(w_old[i_minus:i_plus])
    loc_new[cld(n, 2)] = 0
    wgt_new[cld(n, 2)] = w0

    # positive frequencies
    idx_new = cld(n, 2) + 1
    l_plus = @view l_old[(i_plus + 1):end]
    w_plus = @view w_old[(i_plus + 1):end]
    w_target = sum(w_plus) / (n ÷ 2)
    w = zero(T) # weight
    m = zero(T) # first moment
    for i in eachindex(l_plus)
        # add weight
        w += w_plus[i]
        m += l_plus[i] * w_plus[i]
        if i == lastindex(l_plus)
            # set pole regardless of current weight
            loc_new[idx_new] = m / w
            wgt_new[idx_new] = w
            w = zero(w)
            idx_new == n ||
                throw(ErrorException("failed to discretize positive frequencies"))
        end
        while w > w_target && !(idx_new == n) # outermost pole may overshoot
            δw = w - w_target
            # remove overshoot
            m -= l_plus[i] * δw
            loc_new[idx_new] = m / w_target
            wgt_new[idx_new] = w_target
            idx_new += 1
            # overshoot to next pole
            w = δw
            m = l_plus[i] * δw
            w > 2 * w_target && @warn "degenerate pole at $(l_plus[i])"
        end
    end

    # discretize negative frequencies
    idx_new = n ÷ 2
    l_minus = @view l_old[(i_minus - 1):-1:firstindex(l_old)]
    w_minus = @view w_old[(i_minus - 1):-1:firstindex(w_old)]
    w_target = sum(w_minus) / (n ÷ 2)
    w = zero(T) # weight
    m = zero(T) # first moment
    for i in eachindex(l_minus)
        # add weight
        w += w_minus[i]
        m += l_minus[i] * w_minus[i]
        if i == lastindex(l_minus)
            # set pole regardless of current weight
            loc_new[idx_new] = m / w
            wgt_new[idx_new] = w
            w = zero(w)
            isone(idx_new) ||
                throw(ErrorException("failed to discretize negative frequencies"))
        end
        while w > w_target && !isone(idx_new) # outermost pole may overshoot
            δw = w - w_target
            # remove overshoot
            m -= l_minus[i] * δw
            loc_new[idx_new] = m / w_target
            wgt_new[idx_new] = w_target
            idx_new -= 1
            # overshoot to next pole
            w = δw
            m = l_minus[i] * δw
            w > 2 * w_target && @warn "degenerate pole at $(l_minus[i])"
        end
    end

    return PolesSum(loc_new, wgt_new)
end

"""
    equal_weight_discretization(
        imΔ::AbstractVector{<:Real}, w::AbstractVector{<:Real}, η::Real, n::Int
    )

Discretize the given function `imΔ` on `n` poles, such that each pole has equal weight.

Assumes that `w` is an equidistant grid.
Assumes that `w` has odd number of values.
Assumes that `w` is a symmetric interval.
Assumes that `imΔ` has only semipositive values.
"""
function equal_weight_discretization(
        imΔ::AbstractVector{<:Real}, w::AbstractVector{<:Real}, η::Real, n::Int
    )
    n >= 3 || throw(ArgumentError("at least 3 poles necessary"))
    isodd(n) || throw(ArgumentError("need odd `n`"))

    N = length(w)
    dw = w[2] - w[1]
    M = cld(N, 2)
    total = sum(imΔ) * dw
    v = total / n
    v0 = imΔ[M] * dw
    i = 1
    # pole at w = 0 : weight = min(v, weight in [-η,η])
    while w[M + i] <= η && M + i <= N && v0 < v
        v0 += imΔ[M + i] * dw + imΔ[M - i] * dw
        i += 1
    end
    # calculate total remaining weight for positive frequencies:
    wght_right = sum(imΔ[(M + i):end]) * dw
    # calculate total remaining weight for negative frequencies:
    wght_left = sum(imΔ[1:(M - i)]) * dw
    j = i # remember i for later

    # discretize for positive frequencies:
    P_plus = Float64[]
    V_plus = Float64[]
    vp = 0.0
    pp = 0.0
    v = wght_right / ((n - 1) ÷ 2)
    partial = 0.0
    while partial < wght_right && M + i <= N
        while vp < v && M + i <= N # accumulate weight until v is exceeded or end of grid
            vp += imΔ[M + i] * dw
            pp += w[M + i] * imΔ[M + i] * dw
            i += 1
        end
        δv = vp - v
        if δv > 0 # if vp exceeds v carry over the difference for next pole
            pp -= w[M + i - 1] * δv
            partial += v
            push!(V_plus, v / π)
            push!(P_plus, pp / v)
            vp = δv
            pp = w[M + i - 1] * δv
        elseif vp > 10 * eps()
            # no overshoot, but still some weight
            partial += vp
            push!(V_plus, vp / π)
            push!(P_plus, pp / vp)
            vp = 0.0
            pp = 0.0
        end
    end

    # discretize for negative frequencies:
    P_minus = Float64[]
    V_minus = Float64[]
    vm = 0.0
    pm = 0.0
    v = wght_left / ((n - 1) ÷ 2)
    partial = 0.0
    while partial < wght_left && M - j > 0
        while vm < v && M - j > 0 # accumulate weight until v is exceeded or end of grid
            vm += imΔ[M - j] * dw
            pm += w[M - j] * imΔ[M - j] * dw
            j += 1
        end
        δv = vm - v
        if δv > 0 # if vm exceeds v carry over the difference for next pole
            pm -= w[M - j + 1] * δv
            partial += v
            push!(V_minus, v / π)
            push!(P_minus, pm / v)
            vm = δv
            pm = w[M - j + 1] * δv
        elseif vm > 10 * eps()
            # no overshoot, but still some weight
            partial += vm
            push!(V_minus, vm / π)
            push!(P_minus, pm / vm)
            vm = 0.0
            pm = 0.0
        end
    end
    loc = [reverse!(P_minus); 0; P_plus]
    wgt = [reverse!(V_minus); v0 / π; V_plus]
    return PolesSum(loc, wgt)
end

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
