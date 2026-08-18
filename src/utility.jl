# utility functions

function _with_blas_threads(f, t::Int = 1)
    t_old = BLAS.get_num_threads()
    t_old == t && return f()

    BLAS.set_num_threads(t)
    try
        return f()
    finally
        BLAS.set_num_threads(t_old)
    end
end

"""
    get_RAS_parameters(n_sites::Int, n_occ::Int, n_v_bit::Int, n_c_bit::Int)

Return `n_bit`, `n_v_vector`, `n_c_vector`.
"""
function get_RAS_parameters(n_sites::Int, n_occ::Int, n_v_bit::Int, n_c_bit::Int)
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
    arr = arrowhead_matrix(Δ)
    n_sites = size(arr, 1)
    H_nat, n_occ = to_natural_orbitals(arr)
    n_bit, V_v, V_c = get_RAS_parameters(n_sites, n_occ, L_c, L_v)
    fs = FockSpace(Orbitals(n_bit), FermionicSpin(1 // 2))
    H = natural_orbital_ras_operator(H_nat, H_int, ϵ_imp, fs, n_occ, L_v, L_c, p)
    ψ_start = RASWavefunction_singlet(Dict{UInt64, Float64}, L_v, L_c, V_v, V_c, p)
    E0, ψ0 = ground_state!(H, ψ_start, 5, typemax(Int), var)
    return H, E0, ψ0
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
        H_k::Vector{<:AbstractMatrix},
        Σ_stat::AbstractMatrix,
        Σ_dyn::PolesSumBlock,
        n_fill::Real;
        μ_tol::Real = 1.0e-6,
        b_max::Int = 30,
        μ_min::Real = minimum(locations(Σ_dyn)),
        μ_max::Real = maximum(locations(Σ_dyn)),
        tol_weight::Real = 0,
    )

Find chemical potential ``μ``, such that desired filling ``n_\\mathrm{fill}`` is fulfilled

```math
\\begin{aligned}
n_\\mathrm{fill}
& ≡
∫_{-∞}^0 \\mathrm{d}ω~\\mathrm{Tr}
\\left[
-\\frac{1}{π}\\mathrm{Im}~G_\\mathrm{loc}(ω+\\mathrm{i}0^+)
\\right] \\\\
& =
∫_{-∞}^0 \\mathrm{d}ω~\\mathrm{Tr}
\\left[
-\\frac{1}{π}\\mathrm{Im}~
\\frac{1}{N_k} ∑_k \\frac{1}{ω + \\mathrm{i}0^+ +μ - H_k - Σ(ω + \\mathrm{i}0^+)}
\\right] .
\\end{aligned}
```

A bisection algorithm is used which stops once `Δμ < μ_tol`
or `b_max` iterations are surpassed.

Returns the calculated chemical potential and effective filling.

# Arguments
- `μ_tol::Real = 1.0e-6`: tolerance `Δμ` to exit bisection early
- `b_max::Int = 30`: maximum number of bisections
- `μ_min::Real = minimum(locations(Σ_dyn))`: initial lower bound for `μ`
- `μ_max::Real = maximum(locations(Σ_dyn))`: initial upper bound for `μ`
- `tol_weight::Real = 0`: treat weights less or equal than this value in `Σ_dyn` as zero
"""
function find_chemical_potential(
        H_k::Vector{<:AbstractMatrix},
        Σ_stat::AbstractMatrix,
        Σ_dyn::PolesSumBlock,
        n_fill::Real;
        μ_tol::Real = 1.0e-6,
        b_max::Int = 30,
        μ_min::Real = minimum(locations(Σ_dyn)),
        μ_max::Real = maximum(locations(Σ_dyn)),
        tol_weight::Real = 0,
    )
    # check input
    n_b = size(first(H_k), 1)
    allequal(size, H_k)::Bool || throw(DimensionMismatch("different matrix sizes in H_k"))
    size(Σ_stat) == (n_b, n_b) || throw(DimensionMismatch("size of Σ_stat does not match H_k"))
    (size(Σ_dyn) == (n_b, n_b))::Bool || throw(DimensionMismatch("size of Σ_dyn does not match H_k"))
    μ_min < μ_max || throw(ArgumentError("violating μ_min < μ_max"))

    # represent dynamic part of self-energy as block arrowhead matrix
    Σ_A = arrowhead_matrix(Σ_dyn, sqrt(tol_weight); thin = true)

    # filling for initial guesses
    n_min = _filling_mu(H_k, Σ_stat, Σ_A, μ_min)
    n_max = _filling_mu(H_k, Σ_stat, Σ_A, μ_max)
    n_min <= n_fill <= n_max ||
        throw(ArgumentError("violating n(μ_min) = $(n_min) <= n_fill <= n(μ_max) = $(n_max)"))

    # bisect chemical potential μ
    μ_new = 0.0
    n_new = 0.0
    n_bisect = 0
    for _ in 1:b_max
        n_bisect += 1
        μ_new = 0.5 * (μ_min + μ_max)
        n_new = _filling_mu(H_k, Σ_stat, Σ_A, μ_new)
        n_new > n_fill ? μ_max = μ_new : μ_min = μ_new
        (μ_max - μ_min) < μ_tol && break
    end
    @debug "chemical potential bisection" n_bisect μ_new μ_min μ_max n_fill n_new

    return μ_new, n_new
end

# Calculate filling for given chemical potential μ.
function _filling_mu(H_k, Σ_stat, Σ_A::AbstractMatrix, μ)
    n_b = LinearAlgebra.checksquare(first(H_k)) # number of bands
    z = zero(float(real(eltype(Σ_A))))
    result = Threads.Atomic{typeof(z)}(z)

    Threads.@threads for i in eachindex(H_k)
        foo = copy(Σ_A)
        foo[1:n_b, 1:n_b] = H_k[i]
        foo[1:n_b, 1:n_b] -= μ * I
        foo[1:n_b, 1:n_b] += Σ_stat
        # NOTE: `foo` is a sparse (block arrowhead matrix).
        # One can use Krylov methods to approximate spectrum
        # if full decomposition is too slow.
        F = eigen!(Hermitian(foo))
        n_loc = z # local filling
        @inbounds for j in axes(Σ_A, 2)
            ϵ = F.values[j]
            v = @view F.vectors[1:n_b, j]
            if ϵ < 0
                # Trace of v*v' is sum of values squared.
                n_loc += sum(abs2, v)
            end
        end
        Threads.atomic_add!(result, n_loc)
    end

    return result[] /= length(H_k)
end

function _issorted_and_unique(grid::AbstractVector{<:Real})
    issorted(grid) || throw(ArgumentError("grid is not sorted"))
    allunique(grid) || throw(ArgumentError("grid has degenerate locations"))
    # isequal() treats -0.0 and 0.0 as unequal although both are zero.
    count(iszero, grid) <= 1 || throw(ArgumentError("grid has duplicate zeros"))
    return true
end

"""
    quasiparticle_weight(Σ::PolesSum; tol::Real = 0, λ::Real = 0)

Obtain the quasiparticle weight on the real axis.

```math
\\begin{aligned}
Z &= \\left(1 - \\frac{∂\\mathrm{Re}~Σ(0)}{∂ω}\\right)^{-1} \\\\
  &= \\left(1 - \\sum_i w_i \\frac{∂}{∂ω}\\left.\\frac{1}{ω-a_i}\\right|_{ω=0}\\right)^{-1} \\\\
  &= \\left(1 + \\sum_{w_i ≥ \\mathrm{tol}} w_i\\frac{1}{a_i^2 + λ^2}\\right)^{-1} \\\\
\\end{aligned}
```

Skip all weights ``w_i< \\mathrm{tol}``.
The variable ``λ`` serves as a regularization parameter.
"""
function quasiparticle_weight(Σ::PolesSum; tol::Real = 0, λ::Real = 0)
    tol >= 0 || throw(ArgumentError("tol must be semipositive"))
    λ >= 0 || throw(ArgumentError("λ must be semipositive"))

    deriv = zero(Float64)
    for i in eachindex(Σ)
        weight(Σ, i) < tol && continue
        deriv += weight(Σ, i) / (location(Σ, i)^2 + λ^2)
    end
    return inv(1 + deriv)
end

"""
    quasiparticle_weight_inflections(
        Σ::PolesSum;
        tol::Real = 0,
        λmin::Real = eps(),
        λmax::Real = 1,
    )

Return the inflection points of the regularized quasiparticle weight

```math
Z(λ) = \\left(1 + \\sum_{w_i \\geq \\mathrm{tol}} \\frac{w_i}{a_i^2 + λ^2}\\right)^{-1}
```

in the window ``[λ_\\mathrm{min}, λ_\\mathrm{max}]``.

The inflection points are the solution of ``∂^2Z(λ)/∂λ^2 = 0``,
which are the roots of

```math
(1+S) T + 4 λ^2 T^2 = 4 λ^2 (1+S) U,
```

with the sums

```math
\\begin{aligned}
S &= \\sum_{w_i \\geq \\mathrm{tol}} \\frac{w_i}{a_i^2 + λ^2}, \\\\
T &= \\sum_{w_i \\geq \\mathrm{tol}} \\frac{w_i}{(a_i^2 + λ^2)^2}, \\\\
U &= \\sum_{w_i \\geq \\mathrm{tol}} \\frac{w_i}{(a_i^2 + λ^2)^3}.
\\end{aligned}
```

# Examples
```jldoctest
julia> Σ = PolesSum([1.0], [2.0]);

julia> quasiparticle_weight_inflections(Σ; λmax = 2.0)
1-element Vector{Float64}:
 0.9999999999999998
```

See also [`quasiparticle_weight`](@ref).
"""
function quasiparticle_weight_inflections(
        Σ::PolesSum;
        tol::Real = 0,
        λmin::Real = eps(),
        λmax::Real = 1,
    )

    # check input
    tol >= 0 || throw(ArgumentError("tol must be semipositive"))
    λmax > 0 || throw(ArgumentError("λmax must be positive"))
    λmin > 0 || throw(ArgumentError("λmin must be positive"))
    λmin < λmax || throw(ArgumentError("violating λmin < λmax"))

    TΣ = float(eltype(Σ))

    @inline function derivative_sums(λ)
        S, T, U = zero(TΣ), zero(TΣ), zero(TΣ)
        for (loc, wgt) in Σ
            wgt < tol && continue
            invden = inv(loc^2 + λ^2)
            S += wgt * invden
            T += wgt * invden^2
            U += wgt * invden^3
        end
        return S, T, U
    end

    # residual = 0 are the inflection points
    @inline function residual(λ)
        S, T, U = derivative_sums(λ)
        return (1 + S) * T + 4λ^2 * T^2 - 4λ^2 * (1 + S) * U
    end

    λs = logrange(λmin, λmax; length = 10_000) # enough points per decade

    roots = TΣ[]
    λ_low = λs[1]
    G_low = residual(λ_low)
    @inbounds for λ_high in λs[2:end]
        G_high = residual(λ_high)
        if sign(G_low) != sign(G_high)
            # bisect for higher accuracy
            for _ in 1:100
                λ_mid = (λ_low + λ_high) / 2
                (λ_mid == λ_low || λ_mid == λ_high) && break
                if sign(residual(λ_mid)) == sign(G_low)
                    λ_low = λ_mid
                else
                    λ_high = λ_mid
                end
            end
            push!(roots, (λ_low + λ_high) / 2)
        end
        G_low, λ_low = G_high, λ_high
    end
    return roots
end
