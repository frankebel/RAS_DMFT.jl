"""
    PolesSum{A <: Real, B <: Number} <: AbstractPolesSum{A, B}

Representation of poles on the real axis with locations ``a_i`` of type `A`
and weights ``w_i`` of type `B`

```math
P(z) = ∑_i \\frac{w_i}{z-a_i}.
```

For a block variant see [`PolesSumBlock`](@ref).
"""
struct PolesSum{A <: Real, B <: Number} <: AbstractPolesSum{A, B}
    locations::Vector{A}
    weights::Vector{B}

    function PolesSum{A, B}(locations, weights) where {A, B}
        length(locations) == length(weights) || throw(DimensionMismatch("length mismatch"))
        _issorted_and_unique(locations)
        return new{A, B}(locations, weights)
    end
end

"""
    PolesSum(locs::AbstractVector{A}, wgts::AbstractVector{B}) where {A,B}

Create a new instance of [`PolesSum`](@ref) by supplying the locations `locs`
and weights `wgts`.

```jldoctest
julia> locs = 0:5;

julia> wgts = 5:10;

julia> P = PolesSum(locs, wgts)
PolesSum{Int64, Int64} with 6 poles

julia> locations(P) == locs
true

julia> weights(P) == wgts
true
```
"""
function PolesSum(locs::AbstractVector{A}, wgts::AbstractVector{B}) where {A, B}
    # Check length for permutation below.
    length(locs) == length(wgts) || throw(DimensionMismatch("length mismatch"))

    # sort
    p = sortperm(locs)
    locs = locs[p]
    wgts = wgts[p]

    # Merge degenerate locations.
    loc_out = similar(locs, 0)
    wgt_out = similar(wgts, 0)
    i = 1
    while i <= length(locs)
        l = locs[i]
        w = wgts[i]
        i += 1
        while i <= length(locs) && locs[i] == l
            w += wgts[i]
            i += 1
        end
        push!(loc_out, l)
        push!(wgt_out, w)
    end
    return PolesSum{A, B}(loc_out, wgt_out)
end

PolesSum{A, B}(P::PolesSum) where {A, B} = convert(PolesSum{A, B}, P)

"""
    add_pole_at_zero!(P::PolesSum)

Add a pole at location zero with weight 0 inside `P`.
"""
function add_pole_at_zero!(P::PolesSum)
    if all(!iszero, locations(P))::Bool
        push!(locations(P), 0)
        push!(weights(P), 0)
        sort!(P)
    end
    return P
end

amplitude(P::PolesSum{<:Any, <:Real}, i::Integer) = sqrt(weight(P, i))
function amplitude(P::PolesSum{<:Any, <:Complex}, i::Integer)
    throw(
        DomainError(
            weight(P, i),
            "amplitude for complex weight not defined"
        )
    )
end

function arrowhead_matrix(P::PolesSum)
    amps = amplitudes(P)
    T = promote_type(eltype(P), eltype(amps))
    n = length(P) + 1
    result = zeros(T, n, n)::Matrix{T}

    @inbounds for i in eachindex(P)
        result[i + 1, i + 1] = location(P, i)
    end
    result[1, 2:end] .= amps
    result[2:end, 1] .= amps

    return result
end

function evaluate(P::PolesSum, z::Number)
    result = zero(complex(float(eltype(P))))
    for (ϵ, w) in P
        result += w / (z - ϵ)
    end
    return result
end

function evaluate_gaussian(P::PolesSum, ω::Real, σ::Real)
    result = zero(complex(float(eltype(P))))
    for (ϵ, w) in P
        re = sqrt(2) / σ * dawson((ω - ϵ) / (sqrt(2) * σ))
        im_part = pdf(Normal(ϵ, σ), ω)
        z = re - im * π * im_part
        result += w * z
    end
    return result
end

function filling(P::PolesSum{<:Any, B}, μ::Real = 0) where {B}
    result = zero(Float64) # half weight changes Int → Float

    for (loc, wgt) in P
        if loc < μ
            result += wgt
        elseif loc == μ
            result += 0.5 * wgt
        else
            break
        end
    end
    return result
end

"""
    inverse(P::PolesSum)

Return the Anderson/star decomposition `(a0::Real, D::PolesSum)` by inverting the input:

```math
P(z) = \\frac{1}{z - a0 - D(z)}
```

for a normalized `P`.
"""
function inverse(P::PolesSum)
    isapprox(moment(P, 0), 1; atol = 1000 * eps(float(eltype(P)))) ||
        throw(ArgumentError("P does not have total weight 1"))

    b0, HA = anderson_matrix(P)
    locs = diag(HA)
    a0 = popfirst!(locs)
    wgts = b0 * HA[1, 2:end]
    map!(abs2, wgts)

    return a0, PolesSum(locs, wgts)
end

function merge_degenerate_poles!(P::PolesSum, tol::Real = 0)
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
            wgts[i0] += popat!(wgts, i)
            deleteat!(locs, i)
        end
    end
    # pole(s) at tol → ∞
    i = findfirst(>(0), locs)
    isnothing(i) && (i = lastindex(locs)) # enforce `i` to be a number
    while i < lastindex(locs)
        if locs[i + 1] - locs[i] <= tol
            # merge
            wgts[i] += popat!(wgts, i + 1)
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
            wgts[i - 1] += popat!(wgts, i)
            deleteat!(locs, i - 1) # keep location closer to zero
            i -= 1
        else
            # decrement index
            i -= 1
        end
    end
    return P
end

"""
    merge_negative_weight!(P::PolesSum)

Remove negative weights of `P` such that
the zeroth (total weight) and first moment are conserved.

For each pole with negative weight at index `i`,
the weight is merged to neighbors using the law of levers.

!!! warning
    For negative outermost poles, the first moment is **not conserved**.

Poles with zero weight after merging are automatically removed.
"""
function merge_negative_weight!(P::PolesSum)
    # check input
    m0 = moment(P, 0)
    m0 >= 0 || throw(
        ArgumentError("total weight is negative: $(m0)."),
    )

    locs = locations(P)
    wgts = weights(P)
    i = 1
    while i <= length(P)
        wgts[i] > 0 && (i += 1; continue)

        if iszero(wgts[i])
            deleteat!(locs, i)
            deleteat!(wgts, i)
            continue
        end

        # single pole with negative weight
        i == 1 && length(P) == 1 && throw(
            ArgumentError("Single pole has negative weight $(wgts[1])."),
        )

        # first pole
        if i == 1
            wgts[2] += wgts[1]
            deleteat!(locs, 1)
            deleteat!(wgts, 1)
            continue
        end

        # last pole
        if i == length(P)
            wgts[end - 1] += wgts[end]
            deleteat!(locs, i)
            deleteat!(wgts, i)
            i -= 1
            continue
        end

        # search leftward for compensation
        j = i - 1
        while true
            # lever rule split between j and i+1
            f_left = (locs[i + 1] - locs[i]) / (locs[i + 1] - locs[j])
            f_right = 1 - f_left
            need_left = -f_left * wgts[i]  # amount needed from left

            if wgts[j] >= need_left
                # w_j can fully compensate w_i.
                wgts[j] -= need_left
                wgts[i + 1] += f_right * wgts[i]
                deleteat!(locs, i)
                deleteat!(wgts, i)
                i = j
                break
            else
                # w_j can't fully compensate w_i.
                # Find fraction f ∈ (0, 1) which can be merged such that w_j becomes 0.
                # w_j + f f_l w_i == 0
                wgts[i] += wgts[j] / f_left
                wgts[i + 1] -= f_right / f_left * wgts[j]
                deleteat!(locs, j)
                deleteat!(wgts, j)
                j -= 1
                i -= 1

                # no more left neighbors
                j < 1 && break
            end
        end
    end

    return P
end

function merge_small_weight!(P::PolesSum, tol::Real)
    # check input
    tol >= 0 || throw(ArgumentError("negative tol is invalid"))
    # loop over all poles
    i = 1
    while i <= length(P)
        loc = location(P, i)
        wgt = weight(P, i)
        if wgt > tol
            # enough weight, go to next
            i += 1
            continue
        end
        if i == 1
            # add weight to next pole
            weights(P)[2] += weight(P, 1)
            deleteat!(locations(P), 1)
            deleteat!(weights(P), 1)
        elseif i == length(P)
            # add weight to previous pole
            weights(P)[end - 1] += wgt
            pop!(locations(P))
            pop!(weights(P))
        else
            # split weight such that zeroth and first moment is conserved
            loc_prev = location(P, i - 1)
            loc_next = location(P, i + 1)
            α = (loc_next - loc) / (loc_next - loc_prev)
            weights(P)[i - 1] += α * wgt
            weights(P)[i + 1] += (1 - α) * wgt
            deleteat!(locations(P), i)
            deleteat!(weights(P), i)
        end
    end
    return P
end

function moment(P::PolesSum, n::Int = 0)
    foo = [loc^n * w for (loc, w) in P]
    # sort by abs to guarantee that odd moments are zero for symmetric input
    sort!(foo; by = abs)
    return sum(foo)
end

"""
    spectral_function_loggaussian(P::PolesSum, ω, b::Real)

Calculate the spectral function ``A(ω) = -1/π \\mathrm{Im}[P(ω)]`` with a
lognormal broadening.

Each pole is broadened as in NRG

```math
b_i δ(ω - a_i) → b_i \\frac{\\mathrm{e}^{-b^2/4}}{\\sqrt{π}|a|b}
\\exp\\left(-\\frac{\\ln^2(ω/a_i)}{b^2}\\right).
```

If there is a pole ``a_i = 0``, it is shifted halfway between its neighbors and
each getting half weight

```math
b_i δ(ω) →
  \\frac{b_i}{2} δ\\left(ω - \\frac{a_{i-1}}{2}\\right)
+ \\frac{b_i}{2} δ\\left(ω - \\frac{a_{i+1}}{2}\\right).
```
"""
function spectral_function_loggaussian(P::PolesSum, ω::Real, b::Real)
    result = zero(ω)
    iszero(ω) && return result # no weight at ω == 0
    for (loc, w) in P
        # only contribute weight if ω is on the same side of the real axis
        sign(ω) == sign(loc) || continue
        prefactor = w * exp(-b^2 / 4) / (b * abs(loc) * sqrt(π))
        result += prefactor * exp(-(log(ω / loc) / b)^2)
    end
    return result
end

function spectral_function_loggaussian(P::PolesSum, ω::AbstractVector{<:Real}, b::Real)
    # map for each point in given grid
    return map(i -> spectral_function_loggaussian(P, i, b), ω)
end

weight(P::PolesSum, i::Integer) = weights(P)[i]

function Base.:+(A::PolesSum{LA, WA}, B::PolesSum{LB, WB}) where {LA, WA, LB, WB}
    L = promote_type(LA, LB)
    W = promote_type(WA, WB)
    na, nb = length(A), length(B)
    locs = Vector{L}(undef, na + nb)
    wgts = Vector{W}(undef, na + nb)
    ia = ib = 1
    k = 1
    @inbounds while ia <= na || ib <= nb
        if ia > na
            locs[k] = location(B, ib)
            wgts[k] = weight(B, ib)
            ib += 1
            k += 1
        elseif ib > nb
            locs[k] = location(A, ia)
            wgts[k] = weight(A, ia)
            ia += 1
            k += 1
        elseif location(A, ia) < location(B, ib)
            locs[k] = location(A, ia)
            wgts[k] = weight(A, ia)
            ia += 1
            k += 1
        elseif location(B, ib) < location(A, ia)
            locs[k] = location(B, ib)
            wgts[k] = weight(B, ib)
            ib += 1
            k += 1
        else
            locs[k] = location(A, ia)
            wgts[k] = weight(A, ia) + weight(B, ib)
            ia += 1
            ib += 1
            k += 1
        end
    end
    resize!(locs, k - 1)
    resize!(wgts, k - 1)
    return PolesSum{L, W}(locs, wgts)
end

Base.:-(P::PolesSum{A, B}) where {A, B} = PolesSum{A, B}(copy(locations(P)), -weights(P))

Base.:-(A::PolesSum, B::PolesSum) = +(A, -B)

function Base.convert(::Type{PolesSum{M, N}}, P::PolesSum{A, B}) where {M, N, A, B}
    locs = convert(Vector{M}, locations(P))
    wgts = convert(Vector{N}, weights(P))
    return PolesSum{M, N}(locs, wgts)
end

function Base.copy(P::PolesSum{A, B}) where {A, B}
    return PolesSum{A, B}(copy(locations(P)), copy(weights(P)))
end

function LinearAlgebra.axpby!(α::Number, x::P, β::Number, y::P) where {P <: PolesSum}
    wy = weights(y)
    rmul!(wy, β)

    # add scaled poles of x
    lx = locations(x)
    wx = weights(x)
    ly = locations(y)
    sizehint!(ly, length(y) + length(x))
    sizehint!(wy, length(y) + length(x))
    @inbounds for i in eachindex(x)
        push!(ly, lx[i])
        push!(wy, α * wx[i])
    end
    sort!(y)
    merge_degenerate_poles!(y)

    return y
end
