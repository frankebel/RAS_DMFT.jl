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
            @inline H_A[i, i] = F.values[i - 1]
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
            @inline H_A[i, i] = F.values[i - n_b]
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
    loc, V = LAPACK.stev!('V', Float64.(locations(P)), Float64.(amplitudes(P)))
    amp = scale(P) * view(V, 1, :)
    wgt = map(abs2, amp)
    return PolesSum(loc, wgt)
end

function PolesContinuedFraction(P::PolesSum)
    tol = sqrt(1000 * eps())
    N = length(P)
    A = Diagonal(locations(P)) # Lanczos on this matrix
    # container for Lanzcos
    a = Vector{Float64}(undef, N)
    b = Vector{Float64}(undef, N - 1)
    V = Matrix{Float64}(undef, N, N)
    # normalize starting vector
    v_old = amplitudes(P)
    s = norm(v_old)
    v_old ./= s
    isapprox(s, 1; atol = 100 * eps()) && (s = one(s))
    # first Lanczos step
    V[:, 1] = v_old
    v_new = A * v_old
    a[1] = v_old ⋅ v_new
    axpy!(-a[1], v_old, v_new)

    # successive steps
    for i in 2:N
        b[i - 1] = norm(v_new)
        if b[i - 1] < tol
            # stop early
            @info "Lanczos stopping early: b[$(i - 1)] = $(b[i - 1])."
            deleteat!(a, i:N)
            deleteat!(b, (i - 1):(N - 1))
            break
        end
        rmul!(v_new, inv(b[i - 1]))
        V[:, i] = v_new
        mul!(v_new, A, view(V, :, i)) # new vector
        for j in 1:(i - 1)
            # Gram-Schmidt against all states excluding last
            v_old = view(V, :, j)
            axpy!(-(v_old ⋅ v_new), v_old, v_new)
            axpy!(-(v_old ⋅ v_new), v_old, v_new) # twice as it is unstable
        end
        # Gram-Schmidt against last state
        v_old = view(V, :, i)
        a[i] = v_old ⋅ v_new
        axpy!(-a[i], v_old, v_new)
        axpy!(-(v_old ⋅ v_new), v_old, v_new) # twice as it is unstable
    end
    return PolesContinuedFraction(a, b, s)
end

# block form

function PolesSumBlock(P::PolesContinuedFractionBlock)
    T = eltype(P) <: Real ? Float64 : ComplexF64 # use double precision
    F = eigen!(hermitianpart!(Matrix{T}(Array(P))))
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
    for i in eachindex(P)
        i1 = 1 + (i - 1) * n2
        i2 = i * n2
        Q1[i1:i2, :] = amp[i]
    end
    Q1, scl = _orthonormalize_SVD(Q1)
    _orthonormalize_GramSchmidt!(Q1) # numerical instability
    _orthonormalize_GramSchmidt!(Q1) # numerical instability
    # set small values to zero
    tol = sqrt(eps()) * norm(scl)
    for i in eachindex(scl)
        scl[i] < tol && (scl[i] = 0)
    end
    # block Lanczos with full orthogonalization
    loc, amp, _ = block_lanczos_full_ortho(A, Q1, n1 * n2)
    return PolesContinuedFractionBlock(loc, amp, scl)
end

# block form -> scalar form
"""
    PolesSum(P::PolesSumBlock, i::Integer, j::Integer)

Take the ``P_{i,j}`` element.
"""
function PolesSum(P::PolesSumBlock, i::Integer, j::Integer)
    loc = copy(locations(P))
    wgt = map(m -> m[i, j], weights(P))
    return PolesSum(loc, wgt)
end
