"""
    correlator(
        H::RASOperator, ψ0::RASWavefunction, O::Operator, n_kryl::Int
    )

Calculate the correlator

```math
C(z) = \\left⟨ ψ_0 O^† \\frac{1}{z - H} O ψ_0 \\right⟩.
```
"""
function correlator(
        H::RASOperator, ψ0::RASWavefunction, O::Operator, n_kryl::Int
    )
    v = O * ψ0
    scale = norm(v)
    rmul!(v, inv(scale))
    locations, amplitudes = lanczos(H, v, n_kryl)

    # NOTE: break Lanczos early if amplitude is small
    value, index = findmin(amplitudes)
    @debug "smallest amplitude b=$(value) at index $(index)/$(lastindex(amplitudes))"

    # create block tridiagonal pole representation which is then diagonalized
    PCF = PolesContinuedFraction(locations, amplitudes, scale)
    return PolesSum(PCF)
end

"""
    correlator(
        H::RASOperator, ψ0::RASWavefunction, O::AbstractVector{<:Operator}, n_kryl::Int
    )

Calculate the block correlator

```math
C(z) = \\left⟨ ψ_0 O^† \\frac{1}{z - H} O ψ_0 \\right⟩.
```
"""
function correlator(
        H::RASOperator, ψ0::RASWavefunction, O::AbstractVector{<:Operator}, n_kryl::Int
    )
    V = Matrix{typeof(ψ0)}(undef, 1, length(O)) # 1×n matrix
    for i in eachindex(V)
        V[i] = O[i] * ψ0
    end
    W, scale = _orthonormalize_SVD(V)
    locations, amplitudes = block_lanczos(H, W, n_kryl)

    # create block tridiagonal pole representation which is then diagonalized
    Pt = PolesContinuedFractionBlock(locations, amplitudes, scale)
    return PolesSumBlock(Pt)
end

"""
    correlator_plus(
        H::RASOperator, ψ0::RASWavefunction, O, n_kryl::Int
    )

Calculate the positive spectrum of the (block) correlator.

```math
C^+(z) = \\left⟨ ψ_0 O^† \\frac{1}{z - H} O ψ_0 \\right⟩
```

For a single operator a scalar [`PolesSum`](@ref) is returned,
for a vector of operators a block [`PolesSumBlock`](@ref).

See also [`correlator_minus`](@ref).
"""
function correlator_plus(
        H::RASOperator,
        ψ0::RASWavefunction,
        O::Union{Operator, AbstractVector{<:Operator}},
        n_kryl::Int,
    )
    C = correlator(H, ψ0, O, n_kryl)
    _warn_wrong_sign(C, :plus)

    return C
end

"""
    correlator_minus(
        H::RASOperator, ψ0::RASWavefunction, O, n_kryl::Int
    )

Calculate the negative spectrum of the (block) correlator.

```math
C^-(z) = \\left⟨ ψ_0 O^† \\frac{1}{z + H} O ψ_0 \\right⟩
```

For a single operator a scalar [`PolesSum`](@ref) is returned,
for a vector of operators a block [`PolesSumBlock`](@ref).

See also [`correlator_plus`](@ref).
"""
function correlator_minus(
        H::RASOperator,
        ψ0::RASWavefunction,
        O::Union{Operator, AbstractVector{<:Operator}},
        n_kryl::Int,
    )
    C = correlator(H, ψ0, O, n_kryl)

    map!(-, locations(C)) # flip sign of eigenvalues
    reverse!(C) # order form lowest to highest

    _warn_wrong_sign(C, :minus)
    return C
end

# Warn if the correlator `C` carries spectral weight on poles of the wrong sign.
# Never happens on exact arithmetic.
function _warn_wrong_sign(C::AbstractPolesSum, side::Symbol)
    side in (:plus, :minus) || throw(ArgumentError("`side` must be `:plus` or `:minus`"))
    neg = side === :plus
    locs = locations(C)
    isempty(locs) && return nothing
    f = x -> neg ? x < 0 : x > 0
    f(locs[1]) || f(locs[end]) || return nothing # any wrong-sign pole?
    n = count(f, locs)
    rng = neg ? (1:n) : (lastindex(locs) - n + 1):lastindex(locs)
    weight = sum(tr(weights(C)[i]) for i in rng)
    name = neg ? "C+" : "C-"
    word = neg ? "negative" : "positive"
    @warn "$name has $word spectral weight $weight on $n pole(s)"
    return nothing
end
