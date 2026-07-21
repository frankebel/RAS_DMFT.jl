# Methods used for debugging and development.
# There is almost no reason to ever call this module.

module Debug

using Fermions
using Fermions.Bits
using Fermions.Wavefunctions

export
    # Functions
    colorprint,
    diffkeys,
    excitation!,
    get_excitation,
    mul_excitation

"""
    get_excitation(s::S, m_filled::S, m_empty::S) where {S<:Unsigned}

Return excitation of Slater determinant `s`.

`m_filled`, `m_empty` are the masks of default filled/empty sites
in vector of `CIWavefunction`.
"""
@inline function get_excitation(s::S, m_filled::S, m_empty::S) where {S <: Unsigned}
    return count_ones((s & (m_filled | m_empty) ⊻ m_filled))
end

"""
    excitation!(
        ψ::Wavefunction, m_filled::S, m_empty::S, excitation::Int
    ) where {S<:Unsigned}

In Wavefunction `ψ` remove entries higher than `excitation`.

`m_filled`, `m_empty` are the masks of default filled/empty sites
in vector of `CIWavefunction`.
"""
function excitation!(
        ψ::Wavefunction, m_filled::S, m_empty::S, excitation::Int
    ) where {S <: Unsigned}
    excitation >= 0 || throw(ArgumentError("`excitation` >= 0"))
    for k in keys(ψ)
        get_excitation(k, m_filled, m_empty) <= excitation || delete!(ψ, k)
    end
    return
end

# Colored print of Slater determinant
function colorprint(
        io::IO,
        s::Unsigned,
        nfilled_bit::Int = 0,
        nempty_bit::Int = 0,
        nfilled::Int = 0,
        nempty::Int = 0,
    )
    # test input
    nsites = 2 + nfilled_bit + nempty_bit + nfilled + nempty
    2 * nsites <= bitsize(s) || throw(ArgumentError("insufficient bitsize"))
    nfilled_bit >= 0 || throw(ArgumentError("`nfilled_bit` must be >= 0"))
    nempty_bit >= 0 || throw(ArgumentError("`nempty_bit` must be >= 0"))
    nfilled >= 0 || throw(ArgumentError("`nfilled` must be >= 0"))
    nempty >= 0 || throw(ArgumentError("`nempty` must be >= 0"))

    bs = bitstring(s)
    i = length(bs) - 2 * nsites
    print(io, bs[1:i]) # overflowing bits
    # (n, 2)
    printstyled(io, bs[(i + 1):(i += nempty)]; color = :blue)
    printstyled(io, bs[(i + 1):(i += nfilled)]; color = :cyan)
    printstyled(io, bs[(i + 1):(i += nempty_bit)]; color = :green)
    printstyled(io, bs[(i + 1):(i += nfilled_bit)]; color = :magenta)
    printstyled(io, bs[(i + 1):(i += 2)]; color = :red)
    # (n, 1)
    printstyled(io, bs[(i + 1):(i += nempty)]; color = :blue)
    printstyled(io, bs[(i + 1):(i += nfilled)]; color = :cyan)
    printstyled(io, bs[(i + 1):(i += nempty_bit)]; color = :green)
    printstyled(io, bs[(i + 1):(i += nfilled_bit)]; color = :magenta)
    printstyled(io, bs[(i + 1):(i + 2)]; color = :red)
    return nothing
end

# default stdout
function colorprint(
        s::Unsigned, nfilled_bit::Int = 0, nempty_bit::Int = 0, nfilled::Int = 0, nempty::Int = 0
    )
    return colorprint(stdout, s, nfilled_bit, nempty_bit, nfilled, nempty)
end

"""
    mul_excitation(
        H::Operator, ψ::Wavefunction, m_filled::S, m_empty::S, excitation::Int
    ) where {S}

Calculate `H * ψ` excluding higher excitations than `excitation`.

Assumes that `ψ` only contains Slater determiants with
excitation <= `excitation` in the first place.

`m_filled`, `m_empty` are the masks of default filled/empty sites
in vector of `CIWavefunction`.
"""
function mul_excitation(
        H::Operator, ψ::Wavefunction, m_filled::S, m_empty::S, excitation::Int
    ) where {S <: Unsigned}
    ϕ = Wavefunction(ψ)
    for (k, v) in ψ, t in H.terms
        Fermions.Terms.is_mapped_right(t, k) || continue
        k_new, amp = Fermions.Terms.map_state_right(t, k)
        get_excitation(k_new, m_filled, m_empty) > excitation && continue
        if haskey(ϕ, k_new)
            ϕ[k_new] += amp * v
            iszero(ϕ[k_new]) && delete!(ϕ, k_new)
        else
            ϕ[k_new] = amp * v
        end
    end
    return ϕ
end

# Return all keys not present in other Wavefunction
function diffkeys(
        ϕ1::Wavefunction{<:Any, <:Any, T}, ϕ2::Wavefunction{<:Any, <:Any, T}
    ) where {T}
    length(ϕ1) >= length(ϕ2)
    result = keytype(T)[]
    for k in keys(ϕ1)
        haskey(ϕ2, k) || push!(result, k)
    end
    return Set(result)
end

# # Check sign of amplitude.
# # Check if amplitudes are similar.
# function showdiff(ϕ1::Wavefunction{T}, ϕ2::Wavefunction{T}; kwargs...) where {T}
#     length(ϕ1) == length(ϕ2) || throw(ArgumentError("Wavefunction length mismatch"))
#     for (k, v) in ϕ1
#         try
#             # check anticommutator sign relation
#             sign(v) == sign(ϕ2[k]) || @info "sign mismatch" k ϕ1[k] ϕ2[k]
#             # check if values are similar
#             isapprox(v, ϕ2[k]; kwargs...) || @info "amplitudes" k ϕ1[k] ϕ2[k]
#         catch err
#             isa(err, KeyError) && return KeyError(k)
#         end
#     end
# end

# """
#     highest_amplitudes(ψ::Wavefunction, n_det::Int=20)
#
# Print `n_det` most important amplitudes of `ψ`.
# """
# function highest_amplitudes(ψ::Wavefunction, n_det::Int=20)
#     return display(sort(collect(ψ.terms); by=x -> abs(x[2]), rev=true)[1:n_det])
# end

end
