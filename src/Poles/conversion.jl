# Conversion between different representations.

function anderson_matrix(P::PolesSum)
    N = length(P)

    return _with_blas_threads(Threads.nthreads()) do
        H_LP = Diagonal(locations(P))
        V = amplitudes(P)
        b_0 = norm(V)
        V .*= inv(b_0)

        # Find orthogonal complement to create basis (Eq. A41).
        U1 = [V nullspace(V')]
        h = U1' * H_LP * U1

        # Diagonalize bottom left subspace
        # (1 0   ) * ( h11 h12 ) * ( 1 0  ) = ( h11 h12*U2 )
        # (0 U2' ) * ( h21 h22 )   ( 0 U2 )   ( U2'*h21  d )
        # with d diagonal
        F = eigen!(Hermitian(h[2:end, 2:end]))
        # Anderson matrix
        H_A = zero(h)
        H_A[1, 1] = h[1, 1]
        for i in 2:N
            @inbounds H_A[i, i] = F.values[i - 1]
        end
        h21 = @view h[2:end, 1]
        H_A21 = @view H_A[1, 2:end]
        mul!(H_A21, F.vectors', h21)
        H_A[1, 2:end] .= conj.(H_A21) # H_A12

        b_0, Hermitian(H_A)
    end
end

function anderson_matrix(P::PolesSumBlock)
    N = length(P)
    n_b = size(P, 1) # block size

    return _with_blas_threads(Threads.nthreads()) do
        H_LP = Diagonal(repeat(locations(P); inner = n_b))
        V = vcat(amplitudes(P)...)
        R_a, B_0 = RAS_DMFT._orthonormalize_SVD(V)
        RAS_DMFT._orthonormalize_GramSchmidt!(R_a) # numerical instability
        RAS_DMFT._orthonormalize_GramSchmidt!(R_a) # numerical instability

        # Find orthogonal complement to create basis (Eq. A41).
        U1 = [R_a nullspace(R_a')]
        h = U1' * H_LP * U1

        # Diagonalize bottom left subspace
        # (1 0   ) * ( h11 h12 ) * ( 1 0  ) = ( h11 h12*U2 )
        # (0 U2' ) * ( h21 h22 )   ( 0 U2 )   ( U2'*h21  d )
        # with d diagonal
        F = eigen!(Hermitian(h[(n_b + 1):end, (n_b + 1):end]))
        # Anderson matrix
        H_A = zero(h)
        H_A[1:n_b, 1:n_b] .= @view h[1:n_b, 1:n_b]
        for i in (n_b + 1):(N * n_b)
            @inbounds H_A[i, i] = F.values[i - n_b]
        end
        h21 = @view h[(n_b + 1):end, 1:n_b]
        H_A21 = @view H_A[(n_b + 1):end, 1:n_b]
        mul!(H_A21, F.vectors', h21)
        H_A12 = @view H_A[1:n_b, (n_b + 1):end] # H_A12
        adjoint!(H_A12, H_A21)

        B_0, Hermitian(H_A)
    end
end

# scalar form

function PolesSum(P::PolesContinuedFraction)
    # `SymTridiagonal` often raises `LAPACKException(22)`, therefore call directly
    locs, V = LAPACK.stev!('V', Float64.(locations(P)), Float64.(amplitudes(P)))
    amp = scale(P) * view(V, 1, :)
    wgts = map(abs2, amp)
    return PolesSum(locs, wgts)
end

function PolesContinuedFraction(P::PolesSum)
    # normalize starting vector
    N = length(P)
    s = norm(amplitudes(P))
    v = amplitudes(P) ./ s
    isapprox(s, 1; atol = 100 * eps()) && (s = one(s))
    # Scalar Lanczos is the same as block version with q = 1.
    A, B, _ = block_lanczos_full_ortho(Diagonal(locations(P)), reshape(v, N, 1), N)
    a = map(only, A)
    b = map(only, B)
    return PolesContinuedFraction(a, b, s)
end

# block form

function PolesSumBlock(P::PolesContinuedFractionBlock)
    F = eigen!(Hermitian(tridiagonal_matrix(P)))
    amp = scale(P) * view(F.vectors, 1:size(P, 1), :)
    result = PolesSumBlock(F.values, amp)
    remove_zero_weight!(result)
    return result
end

function PolesContinuedFractionBlock(P::PolesSumBlock)
    n1 = length(P)
    n2 = size(P, 1) # block size
    A = Diagonal(repeat(locations(P); inner = n2)) # block Lanczos with this matrix
    amp = amplitudes(P)
    # Create first block Lanczos vectors.
    # Q1 = vcat(amp...) is type-unstable
    T = eltype(eltype(amp))
    Q1 = Matrix{T}(undef, n1 * n2, n2)
    @inbounds for i in eachindex(P)
        i1 = 1 + (i - 1) * n2
        i2 = i * n2
        Q1[i1:i2, :] = amp[i]
    end
    Q1, scl = _orthonormalize_SVD(Q1)
    _orthonormalize_GramSchmidt!(Q1) # numerical instability
    _orthonormalize_GramSchmidt!(Q1) # numerical instability
    # set small values to zero
    tol = sqrt(eps()) * norm(scl)
    @inbounds for i in eachindex(scl)
        scl[i] < tol && (scl[i] = 0)
    end
    # block Lanczos with full orthogonalization
    locs, amps, _ = block_lanczos_full_ortho(A, Q1, n1 * n2)
    return PolesContinuedFractionBlock(locs, amps, scl)
end

# block form -> scalar form
"""
    PolesSum(P::PolesSumBlock, i::Integer, j::Integer)

Take the ``P_{i,j}`` element.
"""
function PolesSum(P::PolesSumBlock, i::Integer, j::Integer)
    locs = copy(locations(P))
    wgts = map(wgt -> wgt[i, j], weights(P))
    return PolesSum(locs, wgts)
end
