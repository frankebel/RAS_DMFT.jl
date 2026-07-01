# Conversion between different representations.

"""
    anderson_matrix(P::PolesSum)
    anderson_matrix(P::PolesSumBlock)

Calculate scaling factor ``B_0`` and Anderson matrix ``H_\\mathrm{A}``
given a sum of poles.

Reference: https://doi.org/10.48550/arXiv.2605.04974, appendix A3d
"""
function anderson_matrix(P::PolesSum)
    T = eltype(P) <: Real ? Float64 : ComplexF64
    N = length(P)

    H_LP = Diagonal(locations(P))
    V = amplitudes(P)
    b_0 = norm(V)
    V .*= inv(b_0)

    # Find orthogonal complement and create orthonormal basis set
    # using QR decomposition (Eq. A41).
    m = Matrix(one(T) * I(N))
    m[:, 1] = V
    U1 = qr(m).Q
    h = Hermitian(U1' * H_LP * U1)

    # Diagonalize bottom left subspace
    U2 = Matrix(one(T) * I(N))
    F = eigen(Hermitian(h[2:end, 2:end]))
    U2[2:end, 2:end] = F.vectors
    # Anderson matrix
    HA = U2' * h * U2
    HA[2:end, 2:end] = Diagonal(F.values)
    H_A = Hermitian(HA)

    return b_0, H_A
end

function anderson_matrix(P::PolesSumBlock)
    T = eltype(P) <: Real ? Float64 : ComplexF64
    N = length(P)
    n = size(P, 1) # block size

    H_LP = Diagonal(repeat(locations(P); inner = n))
    V::Matrix{T} = vcat(amplitudes(P)...)
    Ra, B_0 = RAS_DMFT._orthonormalize_SVD(V)
    RAS_DMFT._orthonormalize_GramSchmidt!(Ra) # numerical instability
    RAS_DMFT._orthonormalize_GramSchmidt!(Ra) # numerical instability

    # Find orthogonal complement and create orthonormal basis set
    # using QR decomposition (Eq. A41).
    m = Matrix(one(T) * I(N * n))
    m[:, 1:n] = Ra
    U1 = qr(m).Q
    h = Hermitian(U1' * H_LP * U1)

    # Diagonalize bottom left subspace
    U2 = Matrix(one(T) * I(N * n))
    F = eigen(Hermitian(h[(n + 1):end, (n + 1):end]))
    U2[(n + 1):end, (n + 1):end] = F.vectors
    # Anderson matrix
    HA = U2' * h * U2
    HA[(n + 1):end, (n + 1):end] = Diagonal(F.values)
    H_A = Hermitian(HA)

    return B_0, H_A
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
