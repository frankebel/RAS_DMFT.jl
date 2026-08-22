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
    length(locs) == length(wgts) || throw(DimensionMismatch("length mismatch"))
    locs, wgts = _sort_merge_degenerate(locs, wgts)
    return PolesSum{A, B}(locs, wgts)
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
    for (loc, wgt) in P
        result += wgt / (z - loc)
    end
    return result
end

function evaluate_gaussian(P::PolesSum, ω::Real, σ::Real)
    result = zero(complex(float(eltype(P))))
    for (loc, wgt) in P
        result += wgt * _gaussian_broadened(ω, loc, σ)
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
    locs, wgts = _sort_merge_degenerate(
        vcat(locations(A), locations(B)), vcat(weights(A), weights(B))
    )
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
