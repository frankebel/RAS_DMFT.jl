"""
    PolesContinuedFractionBlock{A<:Number,B<:Number} <: AbstractPolesContinuedFraction

Representation of poles on the real axis with as a continued fraction with
locations ``A_i`` of type `A` and amplitudes ``B_i`` of type `B`.

All matrices must be hermitian.
The scale factor ``S`` rescales the whole object.

```math
P(ω) = S \\frac{1}{ω-A_1-B_1\\frac{1}{ω-A_2-…}B_1} S
```
"""
struct PolesContinuedFractionBlock{A <: Number, B <: Number} <: AbstractPolesContinuedFraction
    locations::Vector{Matrix{A}}
    amplitudes::Vector{Matrix{B}}
    scale::Matrix{B}

    function PolesContinuedFractionBlock{A, B}(locations, amplitudes, scale) where {A, B}
        length(locations) == length(amplitudes) + 1 ||
            throw(ArgumentError("length mismatch"))
        # hermitian
        all(ishermitian, locations)::Bool ||
            throw(ArgumentError("locations are not hermitian"))
        all(ishermitian, amplitudes)::Bool ||
            throw(ArgumentError("amplitudes are not hermitian"))
        ishermitian(scale) || throw(ArgumentError("scale is not hermitian"))
        # size
        allequal(size, locations)::Bool ||
            throw(DimensionMismatch("locations do not have matching size"))
        allequal(size, amplitudes)::Bool ||
            throw(DimensionMismatch("amplitudes do not have matching size"))
        size(first(locations)) == size(first(amplitudes)) == size(scale) ||
            throw(DimensionMismatch("matrix size mismatch"))
        return new{A, B}(locations, amplitudes, scale)
    end
end

"""
    PolesContinuedFractionBlock(
        loc::AbstractVector{<:AbstractMatrix{<:A}},
        amp::AbstractVector{<:AbstractMatrix{<:B}},
        [scl::AbstractMatrix{<:B},]
    ) where {A,B}

Create a new instance of [`PolesContinuedFractionBlock`](@ref) by supplying
locations `loc`,
amplitudes `amp`,
and scale `scl`.

By default the scale is set to the identity matrix ``1``.
"""
function PolesContinuedFractionBlock(
        loc::AbstractVector{<:AbstractMatrix{<:A}},
        amp::AbstractVector{<:AbstractMatrix{<:B}},
        scl::AbstractMatrix{<:B},
    ) where {A, B}
    return PolesContinuedFractionBlock{A, B}(loc, amp, scl)
end

# convert type
function PolesContinuedFractionBlock{A, B}(P::PolesContinuedFractionBlock) where {A, B}
    return PolesContinuedFractionBlock{A, B}(
        map(i -> Matrix{A}(i), locations(P)),
        map(i -> Matrix{B}(i), amplitudes(P)),
        Matrix{B}(scale(P)),
    )
end

# scale is identity matrix
function PolesContinuedFractionBlock(
        loc::AbstractVector{<:AbstractMatrix{<:A}}, amp::AbstractVector{<:AbstractMatrix{<:B}}
    ) where {A, B}
    scl = LinearAlgebra.I(size(first(loc), 1))
    return PolesContinuedFractionBlock{A, B}(loc, amp, scl)
end

function evaluate_lorentzian(P::PolesContinuedFractionBlock, ω::Real, δ::Real)
    result = zeros(ComplexF64, size(P))
    loc = Iterators.reverse(locations(P))
    amp = Iterators.reverse(amplitudes(P))
    for (A, B) in zip(loc, amp)
        result = B * inv((ω + im * δ) * I - A - result) * B
    end
    result = scale(P) * inv((ω + im * δ) * I - location(P, 1) - result) * scale(P)
    return result
end

function tridiagonal_matrix(P::PolesContinuedFractionBlock)
    n1 = length(P)
    n2 = size(P, 1)
    n = n1 * n2
    result = zeros(eltype(P), n, n)::Matrix{eltype(P)}
    for i in 1:(n1 - 1)
        i1 = 1 + (i - 1) * n2
        i2 = i * n2
        result[i1:i2, (i1 + n2):(i2 + n2)] = amplitudes(P)[i] # upper diagonal
        result[i1:i2, i1:i2] = location(P, i) # main diagonal
        result[(i1 + n2):(i2 + n2), i1:i2] = amplitudes(P)[i] # lower diagonal
    end
    result[(end - n2 + 1):end, (end - n2 + 1):end] = location(P, n1) # last element
    return result
end

Base.eltype(::Type{<:PolesContinuedFractionBlock{A, B}}) where {A, B} = promote_type(A, B)

Base.size(P::PolesContinuedFractionBlock) = size(scale(P))
Base.size(P::PolesContinuedFractionBlock, i) = size(scale(P), i)
