# utility functions

"""
    init_system(
        Δ::PolesSum,
        H_int::Operator,
        ϵ_imp::Real,
        ϵ_mf::Real,
        L_v::Int,
        L_c::Int,
        p::Int,
        var::Real;
        n_kryl::Int = 50,
        n_max_restart::Int = 20,
    )

Return Hamiltonian, ground state energy, and ground state.

`ϵ_imp` is the bare impurity level entering the Hamiltonian, `ϵ_mf` the mean-field
level defining the basis, see [`natural_impurity_orbital`](@ref).

The ground state is found by [`ground_state!`](@ref)
with `n_kryl` Lanczos vectors per restart and at most `n_max_restart` restarts.
"""
function init_system(
        Δ::PolesSum,
        H_int::Operator,
        ϵ_imp::Real,
        ϵ_mf::Real,
        L_v::Int,
        L_c::Int,
        p::Int,
        var::Real;
        n_kryl::Int = 50,
        n_max_restart::Int = 20,
    )
    H_nat = natural_impurity_orbital(Δ, ϵ_mf)
    fs = FockSpace(Orbitals(2 + L_v + L_c), FermionicSpin(1 // 2))
    H = natural_impurity_orbital_ras_operator(H_nat, H_int, ϵ_imp, fs, L_v, L_c, p)
    ψ_start = RASWavefunction_singlet(
        Dict{UInt64, Float64}, L_v, L_c, H.nfilled, H.nempty, p
    )
    E0, ψ0 = ground_state!(H, ψ_start, n_kryl, n_max_restart, var)
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

`Δ0` is the hybridization at the Fermi level,
i.e. the half-width of the non-interacting impurity resonance,
``Δ_0 = -\\mathrm{Im}~Δ(0 + \\mathrm{i}0^+)``.
"""
function temperature_kondo(U::Real, ϵ::Real, Δ0::Real)
    return sqrt(U * Δ0 / 2) * exp(π * ϵ * (ϵ + U) / (2 * U * Δ0))
end

"""
    find_chemical_potential(
        H_ks::Vector{<:AbstractMatrix},
        Σ_stat::AbstractMatrix,
        Σ_dyn::PolesSumBlock,
        idx::AbstractVector{<:Integer},
        n_fill::Real;
        <keyword arguments>
    )

Find chemical potential ``μ``, such that desired filling ``n_\\mathrm{fill}`` is fulfilled

```math
\\begin{aligned}
n_\\mathrm{fill}
& ≡
∫_{-∞}^0 \\mathrm{d}ω~\\mathrm{Tr}
\\left[
-\\frac{1}{π}\\mathrm{Im}~G_\\mathrm{loc}(ω + \\mathrm{i}0^+)
\\right] \\\\
& =
∫_{-∞}^0 \\mathrm{d}ω~\\mathrm{Tr}
\\left[
-\\frac{1}{π}\\mathrm{Im}~
\\frac{1}{N_k} ∑_k \\frac{1}{ω + \\mathrm{i}0^+ + μ - H_k - Σ_\\mathrm{stat}
- P Σ_\\mathrm{dyn} P^†}
\\right] \\, .
\\end{aligned}
```

The self-energy is split into a static part `Σ_stat`
spanning the whole space,
and dynamic part `Σ_dyn` which is only active on some orbitals `idx`.
Using the projector ``P``,
the expression ``P Σ_\\mathrm{dyn} P^†`` is the upfolded self-energy.

A bisection algorithm is used which stops once `Δμ < μ_tol`
or `b_max` iterations are surpassed.

Returns the calculated chemical potential and effective filling.

# Arguments
- `H_ks::Vector{<:AbstractMatrix}`: Hermitian `n_b × n_b` Hamiltonian per k-point
- `Σ_stat::AbstractMatrix`: static self-energy, `n_b × n_b`
- `Σ_dyn::PolesSumBlock`: dynamic self-energy on the correlated block, `n_c × n_c`
- `idx::AbstractVector{<:Integer}`: the `n_c` correlated orbitals
- `n_fill::Real`: target filling
- `μ_tol::Real = 1.0e-7`: stop once the bisection bracket is narrower than this
- `b_max::Int = 40`: maximum number of bisection steps
- `μ_min::Real = minimum(locations(Σ_dyn))`: lower end of the initial bracket
- `μ_max::Real = maximum(locations(Σ_dyn))`: upper end of the initial bracket
- `n_tol::Real = 1.0e-8`: absolute tolerance of the filling integral
"""
function find_chemical_potential(
        H_ks::Vector{<:AbstractMatrix},
        Σ_stat::AbstractMatrix,
        Σ_dyn::PolesSumBlock,
        idx::AbstractVector{<:Integer},
        n_fill::Real;
        μ_tol::Real = 1.0e-7,
        b_max::Int = 40,
        μ_min::Real = minimum(locations(Σ_dyn)),
        μ_max::Real = maximum(locations(Σ_dyn)),
        n_tol::Real = 1.0e-8,
    )
    # check input
    _check_lattice(H_ks, Σ_stat, Σ_dyn, idx)
    μ_min < μ_max || throw(ArgumentError("violating μ_min < μ_max"))
    n_tol > 0 || throw(DomainError(n_tol, "n_tol is not positive"))

    # μ-independent part of the resolvent
    bands = _projected_eigen(H_ks, Σ_stat, idx)

    # filling for initial guesses
    n_min = _filling_mu(bands, Σ_dyn, μ_min, n_tol)
    n_max = _filling_mu(bands, Σ_dyn, μ_max, n_tol)
    n_min <= n_fill <= n_max || throw(
        ArgumentError(
            lazy"violating n(μ_min) = $(n_min) <= n_fill <= n(μ_max) = $(n_max)",
        ),
    )

    # bisect chemical potential μ
    μ_new = 0.0
    n_new = 0.0
    n_bisect = 0
    for _ in 1:b_max
        n_bisect += 1
        μ_new = 0.5 * (μ_min + μ_max)
        n_new = _filling_mu(bands, Σ_dyn, μ_new, n_tol)
        n_new > n_fill ? μ_max = μ_new : μ_min = μ_new
        (μ_max - μ_min) < μ_tol && break
    end
    @debug "chemical potential bisection" n_bisect μ_new μ_min μ_max n_fill n_new
    (μ_max - μ_min) < μ_tol || @warn(
        "bisection used up `b_max` steps without reaching `μ_tol`",
        b_max, μ_tol, bracket = μ_max - μ_min,
    )

    return μ_new, n_new
end

# Diagonalize `Σ_A` with its top-left `n_b × n_b` block replaced by `H + Σ_stat - μ * I`.
# `Σ_A` has to hold the promoted element type of `H` and `Σ_stat`
function _arrowhead_eigen(
        Σ_A::AbstractMatrix,
        H::AbstractMatrix,
        Σ_stat::AbstractMatrix,
        μ::Real,
        n_b::Int,
    )
    foo = copy(Σ_A)
    # one fused broadcast, no temporaries
    view(foo, 1:n_b, 1:n_b) .= H .+ Σ_stat .- μ * one(H)
    return eigen!(Hermitian(foo))
end

# check lattice input
function _check_lattice(H_ks, Σ_stat, Σ_dyn, idx)
    n_b = LinearAlgebra.checksquare(first(H_ks))
    n_c = length(idx)
    allequal(size, H_ks)::Bool || throw(DimensionMismatch("different matrix sizes in H_ks"))
    size(Σ_stat) == (n_b, n_b) ||
        throw(DimensionMismatch("size of Σ_stat does not match H_ks"))
    (size(Σ_dyn) == (n_c, n_c))::Bool ||
        throw(DimensionMismatch("size of Σ_dyn does not match idx"))
    allunique(idx) || throw(ArgumentError("idx has duplicate orbitals"))
    all(in(1:n_b), idx) || throw(ArgumentError("idx outside range"))
    return n_b, n_c
end

# filling at chemical potential μ using the contour integral of left half-plane
function _filling_mu(bands, Σ_dyn::PolesSumBlock, μ::Real, n_tol::Real)
    Es, Vs = bands
    n_k = length(Es)
    n_b = length(first(Es))
    R = real(eltype(first(Vs)))

    function integrand(y)
        z = im * y
        Σ = evaluate(Σ_dyn, z) # shared by all k-points
        result = Threads.Atomic{R}(zero(R))
        Threads.@threads for i in 1:n_k
            D = Diagonal(inv.(z + μ .- Es[i]))
            DV = D * Vs[i]
            g = Vs[i]' * DV # downfolded G_0
            X = (I - Σ * g) \ (Σ * (Vs[i]' * (D * DV)))
            Threads.atomic_add!(result, real(tr(D) + tr(X)))
        end
        return result[] / n_k
    end

    val, _ = quadgk(integrand, 0, Inf; atol = π * n_tol, maxevals = 10_000)
    return n_b / 2 + val / π
end

function _issorted_and_unique(grid::AbstractVector{<:Real})
    issorted(grid) || throw(ArgumentError("grid is not sorted"))
    allunique(grid) || throw(ArgumentError("grid has degenerate locations"))
    # isequal() treats -0.0 and 0.0 as unequal although both are zero.
    count(iszero, grid) <= 1 || throw(ArgumentError("grid has duplicate zeros"))
    return true
end

# Diagonalize `H_k + Σ_stat` for every k-point and return the eigenvalues `E`
# together with the projection of the eigenvectors to `idx`.
function _projected_eigen(H_ks, Σ_stat, idx)
    T = float(promote_type(eltype(eltype(H_ks)), eltype(Σ_stat)))
    n_k = length(H_ks)
    Es = Vector{Vector{real(T)}}(undef, n_k)
    Vs = Vector{Matrix{T}}(undef, n_k)
    Threads.@threads for i in eachindex(H_ks)
        F = eigen!(Hermitian(Matrix{T}(H_ks[i] + Σ_stat)))
        Es[i] = F.values
        Vs[i] = Matrix(view(F.vectors, idx, :)') # U^† P
    end
    return Es, Vs
end

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
