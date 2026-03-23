"""
    PolesSum{A<:Real,B<:Number} <: AbstractPolesSum

Representation of poles on the real axis with locations ``a_i`` of type `A`
and weights ``w_i`` of type `B`

```math
P(ω) = ∑_i \\frac{w_i}{ω-a_i}.
```

For a block variant see [`PolesSumBlock`](@ref).
"""
struct PolesSum{A <: Real, B <: Number} <: AbstractPolesSum
    locations::Vector{A}
    weights::Vector{B}

    function PolesSum{A, B}(locations, weights) where {A, B}
        length(locations) == length(weights) || throw(DimensionMismatch("length mismatch"))
        return new{A, B}(locations, weights)
    end
end

"""
    PolesSum(loc::AbstractVector{A}, wgt::AbstractVector{B}) where {A,B}

Create a new instance of [`PolesSum`](@ref) by supplying the locations `loc`
and weights `wgt`.

```jldoctest
julia> loc = collect(0:5);

julia> wgt = collect(5:10);

julia> P = PolesSum(loc, wgt)
6-element PolesSum{Int64, Int64}

julia> locations(P) === loc
true

julia> weights(P) === wgt
true
```
"""
function PolesSum(loc::AbstractVector{A}, wgt::AbstractVector{B}) where {A, B}
    return PolesSum{A, B}(loc, wgt)
end

# convert type
function PolesSum{A, B}(P::PolesSum) where {A, B}
    return PolesSum{A, B}(Vector{A}(locations(P)), Vector{B}(weights(P)))
end

"""
    add_pole_at_zero!(P::PolesSum)

Add a pole at location zero with weight 0 inside `P`.
"""
function add_pole_at_zero!(P::PolesSum)
    if all(!iszero, locations(P))
        push!(locations(P), 0)
        push!(weights(P), 0)
        sort!(P)
    end
    return P
end

amplitude(P::PolesSum{<:Any, <:Real}, i::Integer) = sqrt(weight(P, i))

function evaluate_gaussian(P::PolesSum, ω::Real, σ::Real)
    real = zero(ω)
    imag = zero(ω)
    for i in eachindex(P)
        w = weight(P, i)
        real += w * sqrt(2) / (π * σ) * dawson((ω - locations(P)[i]) / (sqrt(2) * σ))
        imag += w * pdf(Normal(locations(P)[i], σ), ω)
    end
    result = real - im * imag
    return π * result # not spectral function
end

function evaluate_lorentzian(P::PolesSum, ω::Real, δ::Real)
    result = zero(complex(ω))
    for i in eachindex(P)
        result += weight(P, i) / (ω + im * δ - locations(P)[i])
    end
    return result
end

function merge_degenerate_poles!(P::PolesSum, tol::Real = 0)
    # check input
    tol >= 0 || throw(ArgumentError("tol must not be negative"))
    issorted(P) || throw(ArgumentError("P must be sorted"))
    # get information from P
    loc = locations(P)
    wgt = weights(P)
    # pole(s) at [-tol, tol]
    idx_zeros = findall(i -> abs(i) <= tol, loc)
    if !isempty(idx_zeros)
        i0 = popfirst!(idx_zeros)
        loc[i0] = 0
        for i in reverse!(idx_zeros)
            wgt[i0] += popat!(wgt, i)
            deleteat!(loc, i)
        end
    end
    # pole(s) at tol → ∞
    i = findfirst(>(0), loc)
    isnothing(i) && (i = lastindex(loc)) # enforce `i` to be a number
    while i < lastindex(loc)
        if loc[i + 1] - loc[i] <= tol
            # merge
            wgt[i] += popat!(wgt, i + 1)
            deleteat!(loc, i + 1) # keep location closer to zero
        else
            # increment index
            i += 1
        end
    end
    # pole(s) at -tol → -∞
    i = findlast(<(0), loc)
    isnothing(i) && (i = firstindex(loc)) # enforce `i` to be a number
    while i > firstindex(loc)
        if loc[i] - loc[i - 1] <= tol
            # merge
            wgt[i - 1] += popat!(wgt, i)
            deleteat!(loc, i - 1) # keep location closer to zero
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

Move negative weights of `P` such that the zeroth moment is conserved
and the first moment changes minimally.
"""
function merge_negative_weight!(P::PolesSum)
    # check input
    issorted(P) || throw(ArgumentError("P is not sorted"))
    allunique(P) || throw(ArgumentError("P has degenerate locations"))
    moment(P, 0) >= 0 || throw(ArgumentError("total weight is negative"))

    loc = locations(P)
    wgt = weights(P)
    for i in eachindex(P)
        weight(P, i) >= 0 && continue # no negative weight, go to next
        if i == length(P)
            # find previous positive weight
            for j in Iterators.reverse(1:(i - 1))
                iszero(wgt[j]) && continue
                if wgt[j] + wgt[end] >= 0
                    # wgt[j] can fully compensate wgt[end]
                    wgt[j] += wgt[end]
                    wgt[end] = 0
                    break
                else
                    # wgt[j] can't fully compensate wgt[end]
                    wgt[end] += wgt[j]
                    wgt[j] = 0
                end
            end
        else
            for j in Iterators.reverse(1:(i - 1))
                # find a previous pole with positive weight
                iszero(wgt[j]) && continue
                # calculate fractions how weight should be split
                f_left = (loc[i + 1] - loc[i]) / (loc[i + 1] - loc[j])
                f_right = 1 - f_left
                if wgt[j] + f_left * wgt[i] >= 0
                    # wgt[j] can fully compensate wgt[i]
                    wgt[j] += f_left * wgt[i]
                    wgt[i + 1] += f_right * wgt[i]
                    wgt[i] = 0
                else
                    # wgt[j] can't fully compensate wgt[i].
                    # Find fraction f ∈ (0, 1) which can be merged such that wgt[j] gets 0 weight.
                    # b_j + f f_l b_i === 0
                    wgt[i] += wgt[j] / f_left
                    wgt[i + 1] -= f_right / f_left * wgt[j]
                    wgt[j] = 0
                end
                if j == 1
                    # no pole with positive weight remaining
                    wgt[i + 1] += wgt[i]
                    wgt[i] = 0
                end
            end
            if wgt[i] <= 0
                # negative weight remaining and no previous weight to compensate
                # move negative weight to next pole
                wgt[i + 1] += wgt[i]
                wgt[i] = 0
            end
        end
    end
    moment(P, 0) >= 0 || throw(ArgumentError("total weight got negative"))
    return P
end

function merge_small_weight!(P::PolesSum, tol::Real)
    # check input
    tol >= 0 || throw(ArgumentError("negative tol is invalid"))
    issorted(P) || throw(ArgumentError("P must be sorted"))
    # loop over all poles
    i = 1
    while i <= length(P)
        loc = locations(P)[i]
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
            loc_prev = locations(P)[i - 1]
            loc_next = locations(P)[i + 1]
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
    foo = map(i -> i[1]^n * i[2], zip(locations(P), weights(P)))
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
    for i in eachindex(P)
        if sign(ω) == sign(locations(P)[i])
            # only contribute weight if on same side of real axis
            w = weight(P, i)
            loc = locations(P)[i]
        else
            # frequency has opposite sign compared to pole location
            continue
        end
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

function Core.Array(P::PolesSum)
    T = eltype(P)
    result = Matrix{T}(Diagonal([0; locations(P)]))
    result[1, 2:end] .= amplitudes(P)
    result[2:end, 1] .= amplitudes(P)
    return result
end

function Base.:+(A::PolesSum, B::PolesSum)
    result = PolesSum([locations(A); locations(B)], [weights(A); weights(B)])
    sort!(result)
    merge_degenerate_poles!(result, 0)
    return result
end

function Base.:-(A::PolesSum, B::PolesSum)
    result = PolesSum([locations(A); locations(B)], [weights(A); -weights(B)])
    sort!(result)
    merge_degenerate_poles!(result, 0)
    return result
end

function Base.copy(P::PolesSum)
    return PolesSum(copy(locations(P)), copy(weights(P)))
end

Base.eltype(::Type{<:PolesSum{A, B}}) where {A, B} = promote_type(A, B)

function Base.inv(P::PolesSum)
    isapprox(moment(P, 0), 1; atol = 1000 * eps()) ||
        throw(ArgumentError("P does not have total weight 1"))

    b0, HA = anderson_matrix(P)
    loc = diag(HA)
    a0 = popfirst!(loc)
    wgt = b0 * HA[1, 2:end]
    map!(i -> abs2(i), wgt)
    return a0, PolesSum(loc, wgt)
end

# create a better show?
Base.show(io::IO, P::PolesSum) = print(io, length(P), "-element ", summary(P))

function LinearAlgebra.axpby!(α::Number, x::P, β::Number, y::P) where {P <: PolesSum}
    wy = weights(y)
    rmul!(wy, β)

    # add scaled poles of x
    lx = locations(x)
    wx = weights(x)
    ly = locations(y)
    sizehint!(ly, length(y) + length(x))
    sizehint!(wy, length(y) + length(x))
    @inline for i in eachindex(x)
        push!(ly, lx[i])
        push!(wy, α * wx[i])
    end
    sort!(y)
    merge_degenerate_poles!(y)

    return y
end

function LinearAlgebra.rmul!(P::PolesSum, α::Number)
    rmul!(weights(P), α)
    return P
end
