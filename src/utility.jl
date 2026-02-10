# utility functions

"""
    get_CI_parameters(n_sites::Int, n_occ::Int, n_v_bit::Int, n_c_bit::Int)

Return `n_bit`, `n_v_vector`, `n_c_vector`.
"""
function get_CI_parameters(n_sites::Int, n_occ::Int, n_v_bit::Int, n_c_bit::Int)
    n_bit = 2 + n_v_bit + n_c_bit
    n_emp = n_sites - n_occ
    n_v_vector = n_occ - 1 - n_v_bit
    n_c_vector = n_emp - 1 - n_c_bit
    return n_bit, n_v_vector, n_c_vector
end

"""
    init_system(
        Δ::PolesSum, H_int::Operator, ϵ_imp::Real, L_v::Int, L_c::Int, p::Int, var::Real
    )

Return Hamiltonian, ground state energy, and ground state.
"""
function init_system(
        Δ::PolesSum, H_int::Operator, ϵ_imp::Real, L_v::Int, L_c::Int, p::Int, var::Real
    )
    arr = Array(Δ)
    n_sites = size(arr, 1)
    H_nat, n_occ = to_natural_orbitals(arr)
    n_bit, V_v, V_c = get_CI_parameters(n_sites, n_occ, L_c, L_v)
    fs = FockSpace(Orbitals(n_bit), FermionicSpin(1 // 2))
    H = natural_orbital_ci_operator(H_nat, H_int, ϵ_imp, fs, n_occ, L_v, L_c, p)
    ψ_start = CIWavefunction_singlet(Dict{UInt64, Float64}, L_v, L_c, V_v, V_c, p)
    E0, ψ0 = ground_state!(H, ψ_start, 5, typemax(Int), var)
    return H, E0, ψ0
end

"""
    δ_gaussian(δ_0::Real, δ_∞::Real, σ::Real, ω::Real)

Return

```math
δ(ω) = δ_∞ + (δ_0 - δ_∞) \\exp\\left(-\\frac{ω^2}{2σ^2}\\right).
```

Allows for variable broadening ``δ``.
"""
function δ_gaussian(δ_0::Real, δ_∞::Real, σ::Real, ω::Real)
    return δ_∞ .+ (δ_0 - δ_∞) * exp(-ω^2 / (2 * σ^2))
end

"""
    temperature_kondo(U::Real, ϵ::Real, Δ0::Real)

Calculate the Kondo temperature for
an interaction `U`,
on-site with energy `ϵ`,
and hybridization `Δ0`.

```math
T_\\mathrm{K} = \\sqrt{\\frac{UΔ_0}{2}} \\exp(\\frac{π ϵ(ϵ+U)}{2UΔ_0})
```
"""
function temperature_kondo(U::Real, ϵ::Real, Δ0::Real)
    return sqrt(U * Δ0 / 2) * exp(π * ϵ * (ϵ + U) / (2 * U * Δ0))
end

"""
    find_chemical_potential(
        W::AbstractVector{<:Number},
        Hk::AbstractVector{<:AbstractMatrix{<:Number}},
        Σ::AbstractVector{<:AbstractMatrix{<:Number}},
        n::Real;
        tol::Real=1e-3, # tolerance Δμ
        n_max::Int=30, # maximum number of steps
    )

Get chemical potential ``μ``, such that desired filling `n` is fulfilled.

```math
∫_{-∞}^0 \\mathrm{d}ω \\mathrm{Tr}\\left[-\\frac{1}{π}\\mathrm{Im}~G(ω)\\right] ≡ n
```

with

```math
G(ω) = \\frac{1}{N} ∑_k [(ω+μ)I - H_k - Σ(ω)]^{-1}
```

If the self-energy `Σ` is included, it is applied to each diagonal entry of `indices`.

A bisection algorithm is used which stops once `Δμ < tol`
or `n_max` iterations are surpassed.

Returns the calculated chemical potential and effective filling.
"""
function find_chemical_potential(
        W::AbstractVector{<:Number},
        Hk::AbstractVector{<:AbstractMatrix{<:Number}},
        Σ::AbstractVector{<:AbstractMatrix{<:Number}},
        n::Real;
        tol::Real = 1.0e-3, # tolerance Δμ
        n_max::Int = 30, # maximum number of steps
    )
    # check input
    nb = LinearAlgebra.checksquare(first(Hk)) # number of bands
    all(i -> size(i) == (nb, nb), Hk) ||
        throw(DimensionMismatch("different matrix sizes in Hk"))
    all(i -> size(i) == (nb, nb), Σ) ||
        throw(DimensionMismatch("different matrix sizes in Σ"))
    length(W) == length(Σ) || throw(ArgumentError("length mismatch: W, Σ"))
    Base.require_one_based_indexing(W)
    Base.require_one_based_indexing(Σ)
    tol > 0 || throw(ArgumentError("negative tolerance tol"))
    n_max > 0 || throw(ArgumentError("negative number of steps n_max"))

    # Loop through each combination of (k, ω) and store the eigenvalues of H_k + Σ(ω).
    ev = Array{ComplexF64}(undef, nb, length(Hk), length(Σ)) # storage for eigenvalues
    Threads.@threads for iΣ in eachindex(Σ)
        # foo = Hk + Σ(ω)
        foo = Matrix{ComplexF64}(undef, nb, nb)
        for iH in eachindex(Hk)
            copyto!(foo, Hk[iH])
            foo .+= Σ[iΣ]
            ev[:, iH, iΣ] .= eigvals!(foo)
        end
    end

    # Use bisection to find μ.
    # For starting values use the quantiles: desired filling ± 10 %
    n_low, n_high = n / nb - 0.1, n / nb + 0.1
    # Avoid filling n ∉ [0, 1]
    n_low < 0 && (n_low = zero(n_low))
    n_high > 1 && (n_high = one(n_low))
    # Get initial values.
    μ_low, μ_high = quantile(vec(real(ev)), (n_low, n_high))
    n_low = _get_filling(W, μ_low, ev)
    n_high = _get_filling(W, μ_high, ev)
    n_low <= n <= n_high || throw(ArgumentError("could not find chemical potential"))
    # Bisection with limited number of steps and tolerance.
    steps = 0
    while steps <= n_max && μ_high - μ_low >= tol
        steps += 1
        μ_new = (μ_low + μ_high) / 2
        n_new = _get_filling(W, μ_new, ev)
        n_new >= n ? μ_high = μ_new : μ_low = μ_new # update
    end
    μ = (μ_high + μ_low) / 2
    filling = _get_filling(W, μ, ev)
    return μ, filling
end

function _get_filling(
        W::AbstractArray{<:Number}, # frequency grid
        μ::Real, # chemical potential
        ev::AbstractArray{<:Number, 3}, # eigenvalues
    )
    # n ∝ ∑_{b,k,ω≤0} (ω+μ-ev_{bkω})^{-1}
    filling = zero(ComplexF64)
    ω0 = findlast(i -> real(i) <= 0, W) # sum all indices ω <= 0
    for iω in 1:ω0
        for k in axes(ev, 2)
            for b in axes(ev, 1)
                @inbounds filling += inv(W[iω] + μ - ev[b, k, iω])
            end
        end
    end
    nk = size(ev, 2) # number of k-points
    dω = real(W[2] - W[1]) # assume equidistant grid
    return -imag(filling) / π / nk * dω
end

"""
    moment(f::AbstractVector{<:Complex}, W::AbstractVector{<:Real}, n::Int=0)

Return the `n`-th moment of `f`

```math
-frac{1}{π} ∫_{-∞}^∞ mathrm{d}ω ω^n mathrm{Im}f(ω).
```

The function `f` is defined over the grid `W`.
"""
function moment(f::AbstractVector{<:Complex}, W::AbstractVector{<:Real}, n::Int = 0)
    # check input
    eachindex(f) == eachindex(W) || throw(ArgumentError("f, W use different indexing"))

    # weight by ω^n
    imf = -imag(f)
    @. imf *= W^n

    # add area by trapezoidal rule
    result = zero(eltype(imf))
    for i in eachindex(imf)
        i == firstindex(imf) && continue # no point to the left
        result += (W[i] - W[i - 1]) * (imf[i - 1] + imf[i]) / 2
    end

    # 1/π factor
    result /= π

    return result
end

function moments(f::AbstractVector{<:Complex}, W::AbstractVector{<:Real}, ns)
    return map(i -> moment(f, W, i), ns)
end

function _derivative(P::PolesSum, ω::R, tol::Real = 0) where {R <: Real}
    tol >= 0 || throw(ArgumentError("tol must be semipositive"))

    result = zero(promote_type(eltype(P), R))
    for i in eachindex(P)
        weight(P, i) < tol && continue
        result -= weight(P, i) / (ω - locations(P)[i])^2
    end
    return result
end

"""
    quasiparticle_weight(Σ::PolesSum, tol::Real=0)

Obtain the quasiparticle weight on the real axis.

```math
Z = \\left(1 - \\frac{∂\\mathrm{Re}~Σ(0)}{∂ω}\\right)^{-1}
```

Skip all weights `< tol`.

See also [`quasiparticle_weight_gaussian`](@ref).
"""
function quasiparticle_weight(Σ::PolesSum, tol::Real = 0)
    tol >= 0 || throw(ArgumentError("tol must be semipositive"))

    foo = _derivative(Σ, zero(eltype(Σ)), tol)
    return inv(1 - foo)
end

"""
    quasiparticle_weight(
        Σ::AbstractVector{Complex}, W::AbstractVector{Real}, n::Integer
    )

Obtain the quasiparticle weight:

```math
Z = \\left(1 - \\frac{\\mathrm{Re}~[Σ(\\mathrm{d}ω)-Σ(-\\mathrm{d}ω)]}{2\\mathrm{d}ω}\\right)^{-1}
```

using the five point midpoint approximation

```math
\\frac{∂\\mathrm{Re}~Σ(0)}{∂ω} ≈ \\mathrm{Re}\\frac{-Σ(ω+2ϵ) + 8Σ(ω+ϵ) - 8Σ(ω-ϵ) + Σ(ω-2ϵ)}{12ϵ},
```

where ``ϵ = W_{i+n} - W{i}`` is the step size of the grid.

Assumes an equidistant grid.
"""
function quasiparticle_weight(
        Σ::AbstractVector{<:Complex}, W::AbstractVector{<:Real}, n::Integer
    )
    issorted(W) || throw(ArgumentError("W must be sorted"))
    length(Σ) == length(W) || throw(ArgumentError("Σ, W must have same length"))
    n >= 0 || throw(ArgumentError("n must be positive"))

    z = searchsortedfirst(W, zero(eltype(W))) # index of ω=0
    dω = W[z + n] - W[z] # assumes equidistant grid

    foo = real(-Σ[z + 2 * n] + 8 * Σ[z + n] - 8 * Σ[z - n] + Σ[z - 2 * n])
    foo /= (12 * dω)
    return inv(1 - foo)
end

"""
    quasiparticle_weight_gaussian(Σ::PolesSum, dω::Real, σ::Real)

Obtain the quasiparticle weight:

```math
Z = \\left(1 - \\frac{\\mathrm{Re}~[Σ(\\mathrm{d}ω)-Σ(-\\mathrm{d}ω)]}{2\\mathrm{d}ω}\\right)^{-1}
```

with Gaussian broadening `σ`.

See also [`quasiparticle_weight`](@ref).
"""
function quasiparticle_weight_gaussian(Σ::PolesSum, dω::Real, σ::Real)
    σ >= 0 || throw(ArgumentError("σ must be positive"))
    dω >= 0 || throw(ArgumentError("dω must be positive"))

    f1 = real(evaluate_gaussian(Σ, -dω, σ))
    f2 = real(evaluate_gaussian(Σ, dω, σ))
    return inv(1 - (f2 - f1) / (2 * dω))
end
