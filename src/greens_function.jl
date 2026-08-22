# Methods related to Green's functions.

# Bethe lattice.

"""
    greens_function_bethe_analytic(z::Number, D::Real=1.0)
    greens_function_bethe_analytic(Z::AbstractVector{<:Number}, D::Real=1.0)

Calculate the Green's function for a Bethe lattice
given a frequency `z` in the upper complex plane,
and half-bandwidth `D`.

```math
G(z) = \\frac{2}{D^2} \\left(z - \\mathrm{sgn}(\\mathrm{Re}(z)) \\sqrt{z^2 - D^2}\\right)
```

with ``\\mathrm{sgn}(0) = \\mathrm{sgn}(0^±)``.
"""
function greens_function_bethe_analytic(z::Number, D::Real = 1.0)
    _check_half_width(D)
    s = (-1)^signbit(real(z)) # sign(0) = sign(0^±)
    return 2 / D^2 * (z - s * sqrt((z + 0.0im)^2 - D^2))
end

function greens_function_bethe_analytic(Z::AbstractVector{<:Number}, D::Real = 1.0)
    return map(z -> greens_function_bethe_analytic(z, D), Z)
end

"""
    greens_function_bethe_simple(n_bath::Int, D::Real=1.0)

Return the [`PolesSum`](@ref) representation of the semicircular density of states
with half-bandwidth `D` on `n_bath` poles.

Poles are found by diagonalizing a tridiagonal matrix with hopping ``t=D/2``.

See also
[`greens_function_bethe_grid`](@ref),
[`greens_function_bethe_equal_weight`](@ref).
"""
function greens_function_bethe_simple(n_bath::Int, D::Real = 1.0)
    # check input
    _check_half_width(D)

    dv = zeros(n_bath)
    ev = fill(D / 2, n_bath - 1) # hopping t = D/2
    H0 = SymTridiagonal(dv, ev)
    a, T = eigen(H0)
    b = abs2.(T[:, 1]) # positive values for simplicity
    return PolesSum(a, b)
end

"""
    greens_function_bethe_grid(grid::AbstractVector{<:Real}, D::Real=1.0)

Return the [`PolesSum`](@ref) representation of the semicircular density of states
with half-bandwidth `D` with poles given in `grid`.

See also
[`greens_function_bethe_simple`](@ref),
[`greens_function_bethe_equal_weight`](@ref).
"""
function greens_function_bethe_grid(grid::AbstractVector{<:Real}, D::Real = 1.0)
    # check input
    _issorted_and_unique(grid)
    s = Semicircle(D)
    locations = Vector(grid)
    weights = _bethe_bisection_weights(locations, x -> cdf(s, x))
    return PolesSum(locations, weights)
end

"""
    greens_function_bethe_grid_hubbard3(
        grid::AbstractVector{<:Real}, U::Real=0.0, D::Real=1.0
    )

Return the [`PolesSum`](@ref) representation of the Hubbard III approximation
with half-bandwidth `D` and poles given in `grid`.

Created using two semicircles at ``±U/2``.
"""
function greens_function_bethe_grid_hubbard3(
        grid::AbstractVector{<:Real}, U::Real = 0.0, D::Real = 1.0
    )
    # check input
    _issorted_and_unique(grid)
    s = Semicircle(D)
    locations = Vector(grid)
    weights = _bethe_bisection_weights(
        locations, x -> cdf(s, x + U / 2) + cdf(s, x - U / 2), 2
    )
    weights ./= 2 # normalize 2 distributions
    return PolesSum(locations, weights)
end

"""
    greens_function_bethe_equal_weight(n_bath::Int, D::Real=1.0)

Return the [`PolesSum`](@ref) representation of the semicircular density of states
with half-bandwidth `D` on `n_bath` poles.

Each pole has the same hybridization ``V^2 = 1/n_b``.

See also
[`greens_function_bethe_simple`](@ref),
[`greens_function_bethe_grid`](@ref).
"""
function greens_function_bethe_equal_weight(n_bath::Int, D::Real = 1.0)
    isodd(n_bath) || throw(DomainError(n_bath, "number of bath sites must be odd"))

    wgt = 1 / n_bath # weight for each pole
    s = Semicircle(D)

    # calculate only negative half, mirror due to symmetry
    q = collect(0:wgt:0.5) # equal weight for each pole
    v = quantile.(Semicircle(D), q) # I_l

    # ϵ_l = 1/wgt ∫_{I_l} dω ω f(ω)
    # trapezoid rule with `n_p` points
    locations = Vector{Float64}(undef, n_bath ÷ 2)
    n_p = 128 # arbitrary number
    for i in eachindex(v)
        i == length(v) && break
        # subtract half of border values
        α = -v[i] * pdf(s, v[i]) - v[i + 1] * pdf(s, v[i + 1])
        α /= 2
        for j in LinRange(v[i], v[i + 1], n_p)
            α += j * pdf(s, j)
        end
        α *= (v[i + 1] - v[i]) / n_p # Δω = I_l/n_p
        locations[i] = α
    end
    locations .*= n_bath # locations .*= 1/wgt

    locations = [locations; 0; -reverse(locations)]
    weights = fill(wgt, n_bath)
    return PolesSum(locations, weights)
end

# Dispersion relation H_k supplied by user.

"""
    greens_function_local(
        H_k::AbstractVector{<:AbstractMatrix},
        μ::Real = 0,
    )

Calculate the non-interacting local Green's function for a dispersion relation ``H_k``.

```math
G_{\\mathrm{loc},0}(z) = \\frac{1}{N_k} ∑_k \\frac{1}{z + μ - H_k}
```

Returns a [`PolesSumBlock`](@ref).
"""
function greens_function_local(
        H_k::AbstractVector{<:AbstractMatrix},
        μ::Real = 0,
    )
    # check input
    isempty(H_k) && throw(ArgumentError("H_k is empty"))
    allequal(size, H_k)::Bool || throw(DimensionMismatch("different matrix sizes in H_k"))
    all(ishermitian, H_k) || throw(ArgumentError("H_k is not Hermitian"))

    T = eltype(eltype(H_k)) <: Real ? Float64 : ComplexF64
    locs = real(T)[]
    wgts = Matrix{T}[]
    for H in H_k
        E, U = eigen(Hermitian(H))
        append!(locs, E)
        for i in axes(U, 2)
            u = view(U, :, i)
            push!(wgts, u * u')
        end
    end

    G_loc = PolesSumBlock(locs, wgts)
    rmul!(G_loc, inv(length(H_k))) # prefactor 1/N_k
    shift_spectrum!(G_loc, μ)
    return G_loc
end

"""
    greens_function_local(
        H_k::Vector{<:AbstractMatrix},
        Σ_stat::AbstractMatrix,
        Σ_dyn::PolesSumBlock,
        μ::Real;
        tol_location::Real = 0,
        tol_weight::Real = 0,
    )

Calculate the interacting local Green's function for a given dispersion relation ``H_k``
and self-energy ``Σ(z)``.

```math
G_\\mathrm{loc}(z) = \\frac{1}{N_k} ∑_k \\frac{1}{z + μ - H_k - Σ(z)} .
```

# Arguments
- `tol_location::Real = 0`: treat locations less or equal than this value in `G_loc` as degenerate
- `tol_weight::Real = 0`: treat weights less or equal than this value in `Σ_dyn` as zero
"""
function greens_function_local(
        H_k::Vector{<:AbstractMatrix},
        Σ_stat::AbstractMatrix,
        Σ_dyn::PolesSumBlock,
        μ::Real;
        tol_location::Real = 0,
        tol_weight::Real = 0,
    )
    # check input
    n_b = LinearAlgebra.checksquare(first(H_k)) # number of bands
    allequal(size, H_k)::Bool || throw(DimensionMismatch("different matrix sizes in H_k"))
    (size(Σ_dyn) == (n_b, n_b))::Bool || throw(DimensionMismatch("matrix size of Σ_dyn does not match H_k"))

    # represent dynamic part of self-energy as block arrowhead matrix
    Σ_A = arrowhead_matrix(Σ_dyn, sqrt(tol_weight); thin = true)

    # Calculate Green's function in pole representation.
    T = float(promote_type(eltype(eltype(H_k)), eltype(Σ_stat), eltype(Σ_dyn)))
    n_k = length(H_k)
    dim = size(Σ_A, 1)
    n_p = n_k * dim # total number of poles
    locs = Vector{real(T)}(undef, n_p)
    amps = Matrix{T}(undef, n_b, n_p)
    locs, amps = let locs = locs, amps = amps, Σ_A = Σ_A, n_b = n_b, dim = dim
        Threads.@threads for i in eachindex(H_k)
            F = _arrowhead_eigen(Σ_A, H_k[i], Σ_stat, μ, n_b)
            idx_low = 1 + dim * (i - 1)
            idx_high = idx_low + dim - 1
            @inbounds locs[idx_low:idx_high] = F.values
            @inbounds amps[:, idx_low:idx_high] = F.vectors[1:n_b, :]
        end
        locs, amps
    end
    G = PolesSumBlock(locs, amps, tol_location)
    rmul!(G, inv(n_k))

    return G
end

function _check_half_width(D::Real)
    D > 0 || throw(DomainError(D, "negative half-bandwidth"))
    return nothing
end

# For each grid point, bisect the interval to its neighbors and compute the pole weight
# as the mass of the distribution with `cdf(Inf) = total`.
function _bethe_bisection_weights(locations, cdf, total = 1)
    n = length(locations)
    weights = similar(locations)
    if n == 1
        weights[1] = total
        return weights
    end
    # wgt ∝ cdf(a_high) - cdf(a_low)
    for i in eachindex(locations)
        if i == 1
            # cdf(-Inf) = 0
            @inbounds weights[i] = cdf(0.5 * (locations[i] + locations[i + 1]))
        elseif i == n
            # cdf(Inf) = total
            @inbounds weights[i] = total - cdf(0.5 * (locations[i - 1] + locations[i]))
        else
            @inbounds weights[i] =
                cdf(0.5 * (locations[i] + locations[i + 1])) -
                cdf(0.5 * (locations[i - 1] + locations[i]))
        end
    end
    return weights
end
