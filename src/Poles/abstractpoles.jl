"""
    AbstractPoles

Supertype which represents a function on the real axis as a collection of poles.
"""
abstract type AbstractPoles end

"""
    amplitude(P::AbstractPoles, i::Integer)

Return the amplitude (`sqrt` of weight) of `P` at index `i`.

See also [`amplitudes`](@ref).
"""
function amplitude end

"""
    amplitudes(P::AbstractPoles)

Return the amplitudes (`sqrt` of weights) of `P`.
"""
function amplitudes end

"""
    evaluate_gaussian(P::AbstractPoles, ω, σ)

Evaluate `P` with Gaussian broadening ``σ``.
"""
function evaluate_gaussian end
function evaluate_gaussian(P::AbstractPoles, ω::AbstractVector{<:Real}, σ)
    # map for each point in given grid
    return map(i -> evaluate_gaussian(P, i, σ), ω)
end

"""
    evaluate_lorentzian(P::AbstractPoles, ω, δ)

Evaluate `P` with Lorentzian broadening ``P(ω + \\mathrm{i}δ)``.
"""
function evaluate_lorentzian end
function evaluate_lorentzian(P::AbstractPoles, ω::AbstractVector{<:Real}, δ)
    # map for each point in given grid
    return map(i -> evaluate_lorentzian(P, i, δ), ω)
end

"""
    location(P::AbstractPoles, i::Integer)

Return the location of `P` at index `i`.

See also [`locations`](@ref).
"""
location(P::AbstractPoles, i::Integer) = locations(P)[i]

"""
    locations(P::AbstractPoles)

Return the locations of `P`.

See also [`location`](@ref).
"""
locations(P::AbstractPoles) = P.locations

"""
    weight(P::AbstractPoles, i::Integer)

Return the weight of `P` at index `i`.

See also [`weights`](@ref).
"""
function weight end

"""
    weights(P::AbstractPoles)

Return the weights of `P`.

See also [`weight`](@ref).
"""
function weights end

Base.isempty(P::AbstractPoles) = iszero(length(P))
Base.length(P::AbstractPoles) = length(locations(P))
