# Manipulation of bit patterns (masks).

"""
    mask_fe(slaterdet::Type{<:Unsigned}, nbit::Int, nfilled::Int, nempty::Int)

Return masks of filled/empty sites in vector of `RASWavefunction`.

# Arguments
- `slaterdet`: Type of the mask.
- `nbit`: Number of sites in bit component of `RASWavefunction`.
- `nfilled`: Number of filled sites in vector of `RASWavefunction`.
- `nempty`: Number of empty sites in vector of `RASWavefunction`.
"""
function mask_fe(slaterdet::Type{<:Unsigned}, nbit::Int, nfilled::Int, nempty::Int)
    # test input
    nsites = nbit + nfilled + nempty
    2 * nsites <= bitsize(slaterdet) || throw(ArgumentError("insufficient bitsize"))
    nbit >= 0 || throw(ArgumentError("`nbit` must be >= 0"))
    nfilled >= 0 || throw(ArgumentError("`nfilled` must be >= 0"))
    nempty >= 0 || throw(ArgumentError("`nempty` must be >= 0"))

    m_filled = least_significant(slaterdet, nfilled) << nbit
    m_filled |= least_significant(slaterdet, nfilled) << (nbit + nsites)
    m_empty = least_significant(slaterdet, nempty) << (nbit + nfilled)
    m_empty |= least_significant(slaterdet, nempty) << (nbit + nfilled + nsites)
    return m_filled, m_empty
end

"""
    slater_start(
        slaterdet::Type{<:Unsigned},
        bi::UInt8,
        nfilled_bit::Int,
        nempty_bit::Int,
        nfilled::Int,
        nempty::Int,
    )

Return Slater determinant with impurity and bath site b filling given by `bi`,
valence bath sites filled fully.

For impurity and b the last four bits of `bi` are taken for occupation.
If `0bwxyz` is given the occupation are taken as:

- `n[b, 2] = w`
- `n[i, 2] = x`
- `n[b, 1] = y`
- `n[i, 1] = z`

E.g. for opposite occupation: `0b1001`, `0b0110`.
"""
function slater_start(
        slaterdet::Type{<:Unsigned},
        bi::UInt8,
        nfilled_bit::Int,
        nempty_bit::Int,
        nfilled::Int,
        nempty::Int,
    )
    # test input
    nbit = 2 + nfilled_bit + nempty_bit
    nvector = nfilled + nempty
    nsites = nbit + nfilled + nempty
    2 * nsites <= bitsize(slaterdet) || throw(ArgumentError("insufficient bitsize"))
    nfilled_bit >= 0 || throw(ArgumentError("`nfilled_bit` must be >= 0"))
    nempty_bit >= 0 || throw(ArgumentError("`nempty_bit` must be >= 0"))
    nfilled >= 0 || throw(ArgumentError("`nfilled` must be >= 0"))
    nempty >= 0 || throw(ArgumentError("`nempty` must be >= 0"))

    b2 = (0b1000 & bi) >> 3
    i2 = (0b0100 & bi) >> 2
    b1 = (0b0010 & bi) >> 1
    i1 = 0b0001 & bi
    result = (b2 == 0 ? zero(slaterdet) : one(slaterdet)) << (nsites + 1)
    result |= (i2 == 0 ? zero(slaterdet) : one(slaterdet)) << nsites
    result |= (b1 == 0 ? zero(slaterdet) : one(slaterdet)) << 1
    result |= (i1 == 0 ? zero(slaterdet) : one(slaterdet))
    m_bit, _ = mask_fe(slaterdet, 2, nfilled_bit, nempty_bit + nvector)
    m_vector, _ = mask_fe(slaterdet, nbit, nfilled, nempty)
    result |= m_vector
    result |= m_bit
    return result
end
