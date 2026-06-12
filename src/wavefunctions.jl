# Methods related to `Fermions.Wavefunctions` module.

"""
    Wavefunction_singlet(
        D::Type, L_v::Integer, L_c::Integer, V_v::Integer, V_c::Integer
    )

Initialize a `Wavefunction` with
singlet in impurity and mirror site
and all valence sites filled.
"""
function Wavefunction_singlet(
        D::Type, L_v::Integer, L_c::Integer, V_v::Integer, V_c::Integer
    )
    K, V = keytype(D), valtype(D)
    s1 = slater_start(K, 0b0110, L_v, L_c, V_v, V_c)
    s2 = slater_start(K, 0b1001, L_v, L_c, V_v, V_c)
    return Wavefunction(Dict(s1 => V(1 / sqrt(2)), s2 => V(1 / sqrt(2))))
end

"""
    CIWavefunction_singlet(
        D::Type, L_v::Integer, L_c::Integer, V_v::Integer, V_c::Integer, excitation::Integer
    )

Initialize a `CIWavefunction` with
singlet in impurity and mirror site
and all valence sites filled.
"""
function CIWavefunction_singlet(
        D::Type, L_v::Integer, L_c::Integer, V_v::Integer, V_c::Integer, excitation::Integer
    )
    K, V = keytype(D), valtype(D)
    s1 = slater_start(K, 0b0110, L_v, L_c, 0, 0)
    s2 = slater_start(K, 0b1001, L_v, L_c, 0, 0)
    v_dim = sum(i -> binomial(2 * (V_v + V_c), i), 0:excitation)
    vec = zeros(V, v_dim)
    vec[1] = 1 / sqrt(2)
    result = CIWavefunction(
        Dict(s1 => copy(vec), s2 => copy(vec)), 2 + L_v + L_c, V_v, V_c, excitation
    )
    return result
end

"""
    ground_state!(
        H::CIOperator,
        ψ_start::CIWavefunction,
        n_kryl::Integer,
        n_max_restart::Integer,
        variance::Real,
    )

Shift ``H → H - E_0`` in-place and return ``E_0``, ``|ψ_0⟩``.

Get approximate ground state and energy using steps of `n_kryl` Krylov cycles
and at most `n_max_restart` restarts.
Calculation is stopped early if `⟨H^2⟩ <= variance`.
"""
function ground_state!(
        H::CIOperator,
        ψ_start::CIWavefunction,
        n_kryl::Integer,
        n_max_restart::Integer,
        variance::Real,
    )
    # check input
    isapprox(norm(ψ_start), 1; atol = 10 * eps()) ||
        throw(ArgumentError("ψ_start is not normalized"))
    n_kryl >= 1 || throw(ArgumentError("n_kryl must be at least 1"))
    n_max_restart >= 1 || throw(ArgumentError("n_max_restart must be at least 1"))
    variance >= 0 || throw(ArgumentError("variance must be >= 0"))

    # initial guess
    ψ0 = copy(ψ_start)
    E0 = dot(ψ0, H, ψ0)
    Fermions.shift_spectrum!(H, E0)

    # containers to reduce allocations
    a = Vector{Float64}(undef, n_kryl)
    b = Vector{Float64}(undef, n_kryl - 1)
    states = [similar(ψ_start) for _ in 1:(n_kryl + 1)]

    for itr in 1:n_max_restart
        lanczos!(a, b, states, H, ψ0, n_kryl)
        F = eigen!(SymTridiagonal(a, b))
        # new state is linear combination
        rmul!(ψ0, F.vectors[1, 1]) # rescale first element
        @inbounds for i in 2:n_kryl
            # add all other elements
            axpy!(F.vectors[i, 1], states[i], ψ0) # ψ0_new += c_i * ψ_i
        end
        normalize!(ψ0) # possible orthogonality loss in Lanczos
        E0 = dot(ψ0, H, ψ0)
        Fermions.shift_spectrum!(H, E0)

        # calculate variance
        foo = H * ψ0
        var = foo ⋅ foo
        if var <= variance
            # variance is below input
            @info "ground state variance reached after $(n_kryl * itr) Krylov steps"
            break
        elseif itr == n_max_restart
            @info "Target variance not reached. Stopped at $(var)"
        end
    end

    # find constant term for E0
    for t in H.opbit.terms
        if iszero(t.mask) &&
                iszero(t.left) &&
                iszero(t.right) &&
                iszero(t.change) &&
                iszero(t.signmask)
            E0 = -t.value
            break
        end
    end

    return E0, ψ0
end
