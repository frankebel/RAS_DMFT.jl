"""
    PolesContinuedFraction{A<:Real,B<:Real} <: AbstractPolesContinuedFraction

Representation of poles on the real axis as a continued fraction with
locations ``a_i`` of type `A` and amplitudes ``b_i`` of type `B`.

The scale factor ``s`` of type `B` rescales the whole object.

```math
P(ω) = \\frac{s^2}{ω - a_1 - \\frac{b_1^2}{ω - a_2 - …}}
```
"""
struct PolesContinuedFraction{A <: Real, B <: Real} <: AbstractPolesContinuedFraction
    locations::Vector{A}
    amplitudes::Vector{B}
    scale::B

    function PolesContinuedFraction{A, B}(locations, amplitudes, scale) where {A, B}
        length(locations) == length(amplitudes) + 1 ||
            throw(ArgumentError("length mismatch"))
        return new{A, B}(locations, amplitudes, scale)
    end
end

"""
    PolesContinuedFraction(
        locs::AbstractVector{A}, amps::AbstractVector{B}, scl=one(B)
    ) where {A,B}

Create a new instance of [`PolesContinuedFraction`](@ref) by supplying locations `locs`,
amplitudes `amps`, and scale `scl`.

By default the scale is set to ``1``.
"""
function PolesContinuedFraction(
        locs::AbstractVector{A}, amps::AbstractVector{B}, scl = one(B)
    ) where {A, B}
    return PolesContinuedFraction{A, B}(locs, amps, scl)
end

# convert type
function PolesContinuedFraction{A, B}(P::PolesContinuedFraction) where {A, B}
    return PolesContinuedFraction{A, B}(
        Vector{A}(locations(P)), Vector{B}(amplitudes(P)), B(scale(P))
    )
end

function evaluate_lorentzian(P::PolesContinuedFraction, ω::Real, δ::Real)
    result = zero(ComplexF64)
    locs = Iterators.reverse(locations(P))
    amps = Iterators.reverse(amplitudes(P))
    for (a, b) in zip(locs, amps)
        result = abs2(b) / (ω + im * δ - a - result)
    end
    result = abs2(scale(P)) / (ω + im * δ - location(P, 1) - result)
    return result
end

tridiagonal_matrix(P::PolesContinuedFraction) = Matrix(SymTridiagonal(P))

weight(P::PolesContinuedFraction, i::Integer) = abs2(amplitudes(P)[i])

weights(P::PolesContinuedFraction) = abs2.(amplitudes(P))

Base.eltype(::Type{<:PolesContinuedFraction{A, B}}) where {A, B} = promote_type(A, B)

function Base.show(io::IO, P::PolesContinuedFraction)
    return print(io, summary(P), " with ", length(P), " poles")
end

function LinearAlgebra.SymTridiagonal(P::PolesContinuedFraction)
    # lose information about scale
    return SymTridiagonal(locations(P), amplitudes(P))
end
