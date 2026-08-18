"""
    block_lanczos(
        H::RASOperator, W::AbstractMatrix{<:WF}, n_kryl::Int
    ) where {WF<:RASWavefunction}


Block Lanczos algorithm for given operator `H`, states in `W` and `n_kryl` Krylov cycles.

Returns main diagonal `A` and subdiagonal `B`.
"""
function block_lanczos(
        H::RASOperator, Q1::AbstractMatrix{<:WF}, N::Integer
    ) where {WF <: RASWavefunction}
    foo, q = size(Q1)
    isone(foo) || throw(ArgumentError("input matrix must have 1 row"))
    A = Vector{Matrix{Float64}}(undef, N)
    B = Vector{Matrix{Float64}}(undef, N - 1)
    N >= 1 || throw(ArgumentError("N must be >= 1"))
    T = scalartype(WF)
    Q_curr = deepcopy(Q1)
    Q_new = similar(Q1)
    Q_old = similar(Q1)
    # containers to reduce allocations
    Q_int = zero(Q1) # int ≙ intermediate
    V1 = Vector{T}(undef, q)
    M1 = Matrix{T}(undef, q, q)

    # first step
    mul!(Q_new, H, Q_curr)
    @inbounds A[1] = Matrix{T}(undef, q, q)
    mul!(A[1], Q_curr', Q_new)
    hermitianpart!(A[1]) # enforce due to finite precision
    copyto!(M1, A[1])
    rmul!(M1, -1)
    mul!(Q_new, Q_curr, M1, true, true)

    # successive steps
    for j in 2:N
        # orthonormalize Q_new
        zerovector!(Q_int)
        @inbounds B[j - 1] = Matrix{T}(undef, q, q)
        _orthonormalize_SVD!(V1, M1, B[j - 1], Q_int, Q_new)
        # TODO: Stop early if norm is small.
        Q_new, Q_int = Q_int, Q_new
        # cycle for new step
        Q_curr, Q_old, Q_new = Q_new, Q_curr, Q_old
        # Q_new = H Q
        zerovector!(Q_new)
        mul!(Q_new, H, Q_curr)
        # Q_new -= Q_old B^†_{j-1}
        adjoint!(M1, B[j - 1])
        rmul!(M1, -1)
        mul!(Q_new, Q_old, M1, true, true)
        # A_j = Q^† Q_new
        @inbounds A[j] = Matrix{T}(undef, q, q)
        mul!(A[j], Q_curr', Q_new)
        hermitianpart!(A[j]) # enforce due to finite precision
        # Q_new -= Q A_j
        copyto!(M1, A[j])
        rmul!(M1, -1)
        mul!(Q_new, Q_curr, M1, true, true)
    end
    return A, B
end

"""
    block_lanczos_full_ortho(
        H::AbstractMatrix{T1}, Q1::Matrix{T2}, N::Integer
    ) where {T1<:Number,T2<:Number}

Block Lanczos ``\\mathrm{span}{Q_1, H Q_1, …, H^N Q_1}`` with full reorthogonalization.

If the off-diagonal element ``B_i`` gets small, the algorithm is stopped early
and not all `N` steps are run.

Returns main diagonal `A`, lower diagonal `B` and states `Q`.
"""
function block_lanczos_full_ortho(
        H::AbstractMatrix{T1}, Q1::Matrix{T2}, N::Integer
    ) where {T1 <: Number, T2 <: Number}
    tol = sqrt(1000 * eps())
    T = promote_type(T1, T2)
    n, q = size(Q1)

    # check input
    ishermitian(H)::Bool || throw(ArgumentError("H is not hermitian"))
    size(H, 2) == n || throw(DimensionMismatch("dimensions of H and Q1 don't match"))
    N >= 1 || throw(ArgumentError("N must be >= 1"))

    A = Vector{Matrix{T}}(undef, N)
    B = Vector{Matrix{T}}(undef, N - 1)
    Q = Vector{Matrix{T}}(undef, N)
    # containers to reduce allocations
    V1 = Vector{real(T)}(undef, q)
    M1 = Matrix{T}(undef, q, q)

    # first step
    Q[1] = copy(Q1)
    Q_new = H * Q1
    A[1] = Q1' * Q_new
    hermitianpart!(A[1])
    mul!(Q_new, Q1, A[1], -1, 1)

    # successive steps
    for j in 2:N
        for k in axes(Q_new, 2)
            # check for vectors with small magnitude and set to zero
            v = view(Q_new, :, k)
            norm(v) < tol && (v .= 0)
        end

        @inbounds B[j - 1] = Matrix{T}(undef, q, q)
        @inbounds Q[j] = Matrix{T}(undef, n, q)
        _orthonormalize_SVD!(V1, M1, B[j - 1], Q[j], Q_new)
        _orthonormalize_GramSchmidt!(Q[j]) # numerical instability
        _orthonormalize_GramSchmidt!(Q[j]) # numerical instability
        if norm(B[j - 1]) < tol
            # stop early
            @info "block Lanczos stopping early: norm(B[$(j - 1)]) = $(norm(B[j - 1]))"
            deleteat!(A, j:N)
            deleteat!(B, (j - 1):(N - 1))
            deleteat!(Q, j:N)
            break
        end
        # new states
        mul!(Q_new, H, Q[j])
        for l in 1:(j - 1)
            # orthogonalize against all previous states excluding last
            # do twice because it is unstable
            _orthogonalize_states!(M1, Q_new, Q[l])
            _orthogonalize_states!(M1, Q_new, Q[l])
        end
        # orthogonalize against last state
        A[j] = Q[j]' * Q_new
        hermitianpart!(A[j])
        mul!(Q_new, Q[j], A[j], -1, 1)
    end

    return A, B, Q
end
