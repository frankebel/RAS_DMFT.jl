# various methods of calculating self-energy from Poles representation

"""
    self_energy_dyson(
        ϵ_imp::Real,
        Δ0::PolesSum,
        G_imp::PolesSum,
        grid::AbstractVector{<:Real}=locations(Δ0),
    )

Calculate the self-energy purely in [`PolesSum`](@ref) representation
using the Dyson equation.

```math
Σ(ω)
= G_{\\mathrm{imp},0}(ω)^{-1} - G_\\mathrm{imp}(ω)^{-1}
= ω - ϵ_\\mathrm{imp} - Δ_0(ω) - G_\\mathrm{imp}(ω)^{-1}
```

Poles with negative weight are moved into neighbors such that the zeroth and first moment
is conserved locally.
"""
function self_energy_dyson(
        ϵ_imp::Real, Δ0::PolesSum, G_imp::PolesSum, grid::AbstractVector{<:Real} = locations(Δ0)
    )

    # invert impurity Green's function
    a0, G_imp_inv = inv(G_imp)

    # Hartree term
    Σ_H = a0 - ϵ_imp

    # sum of poles
    Σ = G_imp_inv - Δ0
    Σ = to_grid(Σ, grid)
    merge_negative_weight!(Σ)
    remove_zero_weight!(Σ)
    return Σ_H, Σ
end

"""
    self_energy_IFG(C::PolesSumBlock, block::Int = 1)

Given a block sum of poles ``C``, calculate the self-energy using the Schur complement
``Σ(ω) = I(ω) - F^\\mathrm{L}(ω) (G(ω))^{-1} F^\\mathrm{R}(ω)``.

The `block` argument chooses either the top left component (default `1`)
or bottom right (`2`) component.

Returns a `PolesSumBlock` object.

Reference: https://doi.org/10.1103/PhysRevB.105.245132
"""
function self_energy_IFG(C::PolesSumBlock, block::Int = 1)
    T = eltype(C) <: Real ? Float64 : ComplexF64
    N = length(C)
    n = size(C, 1) # block size
    iseven(n) || throw(DomainError(n, "block size of C must be even"))
    block in (1, 2) || throw(DomainError(block, "block must be 1 or 2"))

    B0, HA = anderson_matrix(C)

    # decompose scaling matrix
    F = eigen(B0)
    tol = maximum(F.values) * sqrt(eps(real(T)))
    D = similar(F.values)
    map!(λ -> λ >= tol ? 1 / λ : zero(λ), D, F.values)
    B0_inv = Hermitian(F.vectors' * Diagonal(D) * F.vectors) # B0^{-1}
    map!(λ -> λ >= tol ? 1 / λ^2 : zero(λ), D, F.values)
    B0_inv_sqr = Hermitian(F.vectors' * Diagonal(D) * F.vectors) # B0^{-2}

    # extract blocks.
    A1 = Hermitian(HA[1:n, 1:n])
    A = diag(HA)[(n + 1):end]
    B = view(HA, 1:n, (n + 1):(n * N)) # column vectors b_i

    # take inverse
    A1 = B0_inv * A1 * B0_inv
    amp = B0_inv * B
    P = PolesSumBlock(A, amp)

    # take block
    idx = isone(block) ? (1:(n ÷ 2)) : ((n ÷ 2 + 1):n)
    B0_inv_sqr = Hermitian(B0_inv_sqr[idx, idx])
    A1 = Hermitian(A1[idx, idx])
    for i in eachindex(P)
        weights(P)[i] = weight(P, i)[idx, idx]
    end

    # new scaling matrix
    F = eigen(B0_inv_sqr)
    D = Diagonal(similar(F.values))
    map!(λ -> λ >= tol ? 1 / sqrt(λ) : zero(λ), D, F.values)
    B0 = Hermitian(F.vectors' * D * F.vectors) # B0
    A1 = B0 * A1 * B0
    for i in eachindex(P)
        weights(P)[i] = B0 * Hermitian(weight(P, i)) * B0
    end

    # diagonalize
    H = Array(P)
    H[1:(n ÷ 2), 1:(n ÷ 2)] = A1
    F = eigen!(H)
    loc = F.values
    B = view(F.vectors, 1:(n ÷ 2), 1:size(H, 2)) # column vectors b_i
    amp = B0 * B
    P = PolesSumBlock(loc, amp)
    merge_degenerate_poles!(P, zero(real(T)))

    return P
end

"""
    self_energy_IFG(Σ_H::Real, C::PolesSumBlock, W, δ)

Calculate self-energy as
``Σ_z = Σ^\\mathrm{H} + I_z - F^\\mathrm{L}_z (G_z)^{-1} F^\\mathrm{R}_z``
with Lorentzian broadening.

Assumes that `C` is given as

```math
\\begin{pmatrix}
I             & F^\\mathrm{L} \\\\
F^\\mathrm{R} & G
\\end{pmatrix}.
```
"""
function self_energy_IFG_lorentzian(Σ_H::Real, C::PolesSumBlock, W, δ)
    c = evaluate_lorentzian(C, W, δ)
    I = map(c -> c[1, 1], c)
    F_L = map(c -> c[1, 2], c)
    F_R = map(c -> c[2, 1], c)
    G = map(c -> c[2, 2], c)
    return Σ_H .+ I - F_L ./ G .* F_R
end

"""
    self_energy_IFG_gaussian(Σ_H::Real, C::PolesSumBlock, W, σ)

Calculate self-energy as
``Σ_z = Σ^\\mathrm{H} + I_z - F^\\mathrm{L}_z (G_z)^{-1} F^\\mathrm{R}_z``
with Gaussian broadening.

Assumes that `C` is given as

```math
\\begin{pmatrix}
I             & F^\\mathrm{L} \\\\
F^\\mathrm{R} & G
\\end{pmatrix}.
```
"""
function self_energy_IFG_gaussian(Σ_H::Real, C::PolesSumBlock, W, σ)
    c = evaluate_gaussian(C, W, σ)
    I = map(c -> c[1, 1], c)
    F_L = map(c -> c[1, 2], c)
    F_R = map(c -> c[2, 1], c)
    G = map(c -> c[2, 2], c)
    return Σ_H .+ I - F_L ./ G .* F_R
end
