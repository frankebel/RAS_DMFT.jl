"""
    to_natural_orbitals(H::AbstractMatrix, ϵ::Real=1e-8)

Transforms a single particle Hamiltonian `H` to natural orbital basis.

`H[1,1]` is the onsite energy of the impurity.
States with energies `E ∈ (-ϵ, ϵ)` are considered degenerate.
"""
function to_natural_orbitals(H::AbstractMatrix{<:Real}, ϵ::Real = 1.0e-8)
    ishermitian(H) || throw(ArgumentError("`H` not Hermitian"))
    E, T = LAPACK.syev!('V', 'U', copy(H))
    n_lower = count(<=(-ϵ), E)
    n_zero = count(x -> abs(x) < ϵ, E)
    isodd(n_zero) && @warn "odd number of energies equal zero"
    n_occ = n_lower + (n_zero ÷ 2) # half of states around zero occupied
    if n_zero > 1
        # symmetrize degenerate states around zero
        @info "degenerate zero-energy"
        p = T'[n_occ:(n_occ + 1), 1]
        p ./= norm(p)
        R = inv([p[1] p[2]; p[2] -p[1]]) * [1 / sqrt(2), 1 / sqrt(2)]
        t1 = R[1] * T[:, n_occ] + R[2] * T[:, n_occ + 1]
        t2 = R[1] * T[:, n_occ + 1] - R[2] * T[:, n_occ]
        T[:, n_occ] .= t1[:]
        T[:, n_occ + 1] .= t2[:]
    end

    base_occ = T[:, 1:n_occ] # valence states
    base_emp = T[:, (n_occ + 1):end] # conduction states

    P = base_occ * base_occ' # projector on valence
    Q = I - P # projector on conduction
    w = P[:, 1] # impurity projected on valence
    α = norm(w)
    w ./= α
    u = Q[:, 1] # impurity projected on conduction
    β = norm(u)
    u ./= β

    # orthogonalize the remaining valence wrt impurity
    for v in eachcol(base_occ)
        v .-= dot(v, w) .* w
    end
    # orthogonalize the remaining conduction wrt impurity
    for v in eachcol(base_emp)
        v .-= dot(v, u) .* u
    end

    # Löwdin on valence states
    base_occ[:, 1] .= w
    B = @view base_occ[:, 2:end]
    S = B' * B
    E, V = LAPACK.syev!('V', 'U', S)
    tol_v = maximum(E) * sqrt(eps(typeof(ϵ))) # set all eigenvalues smaller to zero
    R = V * Diagonal(map(x -> x >= tol_v ? 1 / sqrt(x) : 0, E)) * V'
    B .= B * R
    h = base_occ' * H * base_occ

    a_occ, b_occ, _ = sytrd!('L', h)

    # Löwdin on conduction states
    base_emp[:, 1] .= u
    B = @view base_emp[:, 2:end]
    S = B' * B
    E, V = LAPACK.syev!('V', 'U', S)
    tol_c = maximum(E) * sqrt(eps(typeof(ϵ))) # set all eigenvalues smaller to zero
    R = V * Diagonal(map(x -> x >= tol_c ? 1 / sqrt(x) : 0, E)) * V'
    B .= B * R
    h = base_emp' * H * base_emp

    a_emp, b_emp, _ = sytrd!('L', h)

    push!(b_occ, 0.0)
    a = vcat(a_occ, a_emp)
    b = vcat(b_occ, b_emp)
    H_tri = SymTridiagonal(a, b)
    n = length(a)
    v1 = zeros(n)
    v1[1] = α
    v1[n_occ + 1] = β
    v2 = zeros(n)
    v2[1] = β
    v2[n_occ + 1] = -α
    T = Matrix{Float64}(I, n, n)
    T[:, 1] = v1
    T[:, n_occ + 1] = v2
    H_trafo = T' * H_tri * T
    H_trafo = 0.5 * (H_trafo' + H_trafo) # take the Hermitian part
    return H_trafo, n_occ
end

"""
    natural_orbital_operator(
        H_nat::Matrix{T},
        H_int::Operator,
        ϵ_imp::T,
        fock_space::FockSpace,
        n_occ::Int,
        n_v_bit::Int=1,
        n_c_bit::Int=1,
    ) where {T<:Real}

Convert natural orbital Hamiltonian `H_nat` to `Operator`.

# Arguments
- `H_nat::Matrix{T}`: natural orbital Hamiltonian
- `H_int::Operator`: interacting Hamiltonian
- `U::T`: Coulomb repulsion on impurity
- `ϵ_imp::T`: on-site energy of impurity
- `fock_space::FockSpace`: Fock Space used for the system
- `n_occ::Int`: number of occupied sites
- `n_v_bit::Int=1`: number of valence bath sites in bit component
- `n_c_bit::Int=1`: number of conduction bath sites in bit component
"""
function natural_orbital_operator(
        H_nat::Matrix{T},
        H_int::Operator,
        ϵ_imp::T,
        fock_space::FockSpace,
        n_occ::Int,
        n_v_bit::Int = 1,
        n_c_bit::Int = 1,
    ) where {T <: Real}
    ishermitian(H_nat) || throw(ArgumentError("H_nat not Hermitian"))
    n = size(H_nat, 1)
    n_emp = n - n_occ
    n_v = n_occ - 1
    n_c = n_emp - 1
    0 < n_v_bit <= n_v || throw(ArgumentError(lazy"violating 0 < $(n_v_bit) <= $(n_v)"))
    0 < n_c_bit <= n_c || throw(ArgumentError(lazy"violating 0 < $(n_c_bit) <= $(n_c)"))
    c = annihilators(fock_space)
    # impurity i and mirror bath site b
    H = _add_impurity_terms(H_int, c, H_nat, ϵ_imp, n_occ)
    for σ in axes(c, 2)
        # valence bath sites
        for i in 1:n_v
            j = 1 + i # Index in H_nat.
            # Let `k` be the index in the bit component.
            if i <= n_v_bit
                k = 2 + i
            else
                k = 2 + n_c_bit + i
            end
            # bath site
            H += H_nat[j, j] * c[k, σ]' * c[k, σ]
            if i == 1
                # hopping v_1 <-> i
                H += H_nat[j, 1] * c[1, σ]' * c[k, σ]
                H += H_nat[j, 1] * c[k, σ]' * c[1, σ]
                # hopping v_1 <-> b
                H += H_nat[j, n_occ + 1] * c[2, σ]' * c[k, σ]
                H += H_nat[j, n_occ + 1] * c[k, σ]' * c[2, σ]
            elseif i == n_v_bit + 1
                # hopping v_(n_v_bit + 1) <-> v_(n_v_bit)
                H += H_nat[j, j - 1] * c[k, σ]' * c[2 + n_v_bit, σ]
                H += H_nat[j, j - 1] * c[2 + n_v_bit, σ]' * c[k, σ]
            else
                # hopping to previous neighbor
                H += H_nat[j - 1, j] * c[k - 1, σ]' * c[k, σ]
                H += H_nat[j - 1, j] * c[k, σ]' * c[k - 1, σ]
            end
        end
        # conduction bath sites
        for i in 1:n_c
            j = n_occ + 1 + i # Index in H_nat.
            # Let k be the index in the bit component.
            if i <= n_c_bit
                k = 2 + n_v_bit + i
            else
                k = j
            end
            # bath site
            H += H_nat[j, j] * c[k, σ]' * c[k, σ]
            if i == 1
                # hopping c_1 <-> i
                H += H_nat[j, 1] * c[1, σ]' * c[k, σ]
                H += H_nat[j, 1] * c[k, σ]' * c[1, σ]
                # hopping c_1 <-> b
                H += H_nat[j, n_occ + 1] * c[2, σ]' * c[k, σ]
                H += H_nat[j, n_occ + 1] * c[k, σ]' * c[2, σ]
            elseif i == n_c_bit + 1
                # hopping c_(n_c_bit + 1) <-> c_(n_c_bit)
                H += H_nat[j, j - 1] * c[k, σ]' * c[2 + n_v_bit + n_c_bit, σ]
                H += H_nat[j, j - 1] * c[2 + n_v_bit + n_c_bit, σ]' * c[k, σ]
            else
                # hopping to previous neighbor
                H += H_nat[j - 1, j] * c[k - 1, σ]' * c[k, σ]
                H += H_nat[j - 1, j] * c[k, σ]' * c[k - 1, σ]
            end
        end
    end
    return H
end

"""
    natural_orbital_ras_operator(
        H_nat::Matrix{T},
        H_int::Operator,
        ϵ_imp::T,
        fock_space::FockSpace,
        n_occ::Int,
        n_v_bit::Int=1,
        n_c_bit::Int=1,
        excitation::Int=1,
    ) where {T<:Real}

Convert natural orbital Hamiltonian `H_nat` to `RASOperator`.

# Arguments
- `H_nat::Matrix{T}`: natural orbital Hamiltonian
- `H_int::Operator`: interacting Hamiltonian
- `ϵ_imp::T`: on-site energy of impurity
- `fock_space::FockSpace`: Fock Space used for the system
- `n_occ::Int`: number of occupied sites
- `n_v_bit::Int=1`: number of valence bath sites in bit component
- `n_c_bit::Int=1`: number of conduction bath sites in bit component
- `excitation::Int=1`: maximum excitation in bit component

See also `RASOperator`.
"""
function natural_orbital_ras_operator(
        H_nat::Matrix{T},
        H_int::Operator,
        ϵ_imp::T,
        fock_space::FockSpace,
        n_occ::Int,
        n_v_bit::Int = 1,
        n_c_bit::Int = 1,
        excitation::Int = 1,
    ) where {T <: Real}
    # Check if function for zero chain length should be used.
    n_v_bit === n_c_bit === 0 && return _natural_orbital_ras_operator_zero(
        H_nat, H_int, ϵ_imp, fock_space, n_occ, excitation
    )
    # check input
    ishermitian(H_nat) || throw(ArgumentError("H_nat not Hermitian"))
    nflavours(fock_space) >= 2 + n_v_bit + n_c_bit ||
        throw(ArgumentError("fock_space too small"))
    n_v_bit >= 1 || throw(ArgumentError("invalid n_v_bit"))
    n_c_bit >= 1 || throw(ArgumentError("invalid n_c_bit"))
    excitation >= 0 || throw(ArgumentError("negative excitation"))
    n = size(H_nat, 1)
    n_emp = n - n_occ
    n_bit = 2 + n_v_bit + n_c_bit
    n_v = n_occ - 1
    n_c = n_emp - 1
    n_v_bit < n_v || throw(ArgumentError("n_v_bit too big"))
    n_c_bit < n_c || throw(ArgumentError("n_c_bit too big"))
    c = annihilators(fock_space)

    # Create Bitoperator H_bit.
    H_bit = _add_impurity_terms(H_int, c, H_nat, ϵ_imp, n_occ)
    for σ in axes(c, 2)
        # valence bath sites
        H_bit = _add_bath_chain(H_bit, c, σ, H_nat, n_occ, 1, 2, n_v_bit)
        # conduction bath sites
        H_bit = _add_bath_chain(H_bit, c, σ, H_nat, n_occ, n_occ + 1, 2 + n_v_bit, n_c_bit)
    end

    # Create VectorOperator
    n_v_vector = n_v - n_v_bit
    n_c_vector = n_c - n_c_bit
    foo = diag(H_nat)
    esite = [foo[(2 + n_v_bit):n_occ]; foo[(n_occ + n_c_bit + 2):end]]
    foo = diag(H_nat, 1)
    ehop = [foo[(2 + n_v_bit):n_occ]; foo[(n_occ + n_c_bit + 2):end]]

    # Create MixedOperator
    # (i, j, amp)
    mixed = (
        # valence bath site
        (2 + n_v_bit, 1, foo[1 + n_v_bit]),
        # conduction bath site
        (2 + n_v_bit + n_c_bit, n_v_vector + 1, foo[n_occ + n_c_bit + 1]),
    )

    return RASOperator(H_bit, mixed, esite, ehop, n_bit, n_v_vector, n_c_vector, excitation)
end

# Same as `natural_orbital_operator` but with `n_v_bit === n_c_bit === 0`.
function _natural_orbital_ras_operator_zero(
        H_nat::Matrix{T},
        H_int::Operator,
        ϵ_imp::T,
        fock_space::FockSpace,
        n_occ::Int,
        excitation::Int = 1,
    ) where {T <: Real}
    ishermitian(H_nat) || throw(ArgumentError("H_nat not Hermitian"))
    nflavours(fock_space) >= 2 || throw(ArgumentError("fock_space too small"))
    excitation >= 0 || throw(ArgumentError("negative excitation"))
    n = size(H_nat, 1)
    n_emp = n - n_occ
    n_bit = 2
    n_v = n_occ - 1
    n_c = n_emp - 1
    c = annihilators(fock_space)

    # Create Bitoperator H_bit.
    H_bit = _add_impurity_terms(H_int, c, H_nat, ϵ_imp, n_occ)

    # Create VectorOperator
    n_v_vector = n_v
    n_c_vector = n_c
    foo = diag(H_nat)
    esite = [foo[2:n_occ]; foo[(n_occ + 2):end]]
    foo = diag(H_nat, 1)
    ehop = [foo[2:n_occ]; foo[(n_occ + 2):end]]

    # Create MixedOperator
    # (i, j, amp)
    mixed = (
        (1, 1, H_nat[1, 2]), # i <-> v1
        (2, 1, H_nat[n_occ + 1, 2]), # b <-> v1
        (1, n_v_vector + 1, H_nat[1, n_occ + 2]), # i <-> c1
        (2, n_v_vector + 1, H_nat[n_occ + 1, n_occ + 2]), # b <-> c1
    )

    return RASOperator(H_bit, mixed, esite, ehop, n_bit, n_v_vector, n_c_vector, excitation)
end

# Add the impurity on-site energy and the impurity–mirror-site hopping.
# This block is identical for the full and RAS operator constructions.
function _add_impurity_terms(H, c, H_nat::AbstractMatrix, ϵ_imp, n_occ::Int)
    for σ in axes(c, 2)
        # impurity i
        H += ϵ_imp * c[1, σ]' * c[1, σ]
        # mirror bath site b
        H += H_nat[n_occ + 1, n_occ + 1] * c[2, σ]' * c[2, σ]
        # hopping i <-> b
        H += H_nat[n_occ + 1, 1] * c[1, σ]' * c[2, σ]
        H += H_nat[n_occ + 1, 1] * c[2, σ]' * c[1, σ]
    end
    return H
end

# Add the on-site energies and nearest-neighbor hopping of one bath-site chain
# (valence or conduction) inside the RAS bit component.
function _add_bath_chain(H, c, σ, H_nat::AbstractMatrix, n_occ::Int, base_nat::Int, base_bit::Int, n_chain::Int)
    for i in 1:n_chain
        site_nat = base_nat + i # site index in H_nat
        site_bit = base_bit + i # site index in bit component
        # bath site
        H += H_nat[site_nat, site_nat] * c[site_bit, σ]' * c[site_bit, σ]
        if i == 1
            # hopping chain start <-> i and chain start <-> b
            H += H_nat[site_nat, 1] * c[1, σ]' * c[site_bit, σ]
            H += H_nat[site_nat, 1] * c[site_bit, σ]' * c[1, σ]
            H += H_nat[site_nat, n_occ + 1] * c[2, σ]' * c[site_bit, σ]
            H += H_nat[site_nat, n_occ + 1] * c[site_bit, σ]' * c[2, σ]
        else
            # hopping to previous neighbor
            H += H_nat[site_nat - 1, site_nat] * c[site_bit - 1, σ]' * c[site_bit, σ]
            H += H_nat[site_nat - 1, site_nat] * c[site_bit, σ]' * c[site_bit - 1, σ]
        end
    end
    return H
end
