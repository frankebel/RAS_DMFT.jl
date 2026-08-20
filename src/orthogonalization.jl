# Methods for orthogonalization and orthonormalization of states.

# Explanation available under ch. 3.2 of Martin's thesis.
# https://doi.org/10.11588/heidok.00029305
function _orthonormalize_SVD!(
        # user supplies all containers to calculate in-place
        V1::AbstractVector{<:Real}, # container for Λ^{±1/2}
        M1::AbstractMatrix{<:T}, # container
        S_sqrt::AbstractMatrix{<:T}, # store S^{1/2}
        Q_new::AbstractMatrix, # store orthonormal states
        Q::AbstractMatrix, # states to orthonormalize
    ) where {T <: Number}
    mul!(M1, Q', Q) # overlap matrix
    F = eigen(hermitianpart!(M1))
    tol = maximum(F.values) * sqrt(eps(real(T)))
    # orthonormalize states
    map!(λ -> λ >= tol ? 1 / sqrt(λ) : zero(λ), V1, F.values) # Λ^{-1/2}
    mul!(M1, Diagonal(V1), F.vectors')
    mul!(S_sqrt, F.vectors, M1)
    hermitianpart!(S_sqrt) # S^{-1/2}
    mul!(Q_new, Q, S_sqrt) # Q_new = Q S^{-1/2}
    # B = S^{1/2}
    map!(λ -> λ >= tol ? sqrt(λ) : zero(λ), V1, F.values) # Λ^{1/2}
    mul!(M1, Diagonal(V1), F.vectors')
    mul!(S_sqrt, F.vectors, M1) # S^{1/2}
    hermitianpart!(S_sqrt)
    return nothing
end

"""
    _orthonormalize_SVD(Q::AbstractMatrix)

Löwdin orthonormalization (Singular value decomposition SVD) for given states `Q`.

Calculate overlap matrix ``S = Q^† Q`` and diagonalize

```math
\\begin{aligned}
S        &= U Λ U^† \\\\
S^{1/2}  &= U Λ^{1/2} U^† \\\\
S^{-1/2} &= U Λ^{-1/2} U^†.
\\end{aligned}
```

Objects of interest are ``Q S^{-1/2}`` and ``S^{1/2}``.
"""
function _orthonormalize_SVD(Q::AbstractMatrix)
    q = size(Q, 2)
    T = _scalartype(eltype(Q))
    Q_new = similar(Q)
    V1 = Vector{real(T)}(undef, q)
    M1 = Matrix{T}(undef, q, q)
    S_sqrt = similar(M1)
    _orthonormalize_SVD!(V1, M1, S_sqrt, Q_new, Q)
    return Q_new, S_sqrt
end

_scalartype(::Type{T}) where {T <: Number} = T
_scalartype(::Type{WF}) where {WF <: RASWavefunction} = scalartype(WF)

"""
    orthonormalize_GramSchmidt!(V::AbstractMatrix{<:Number})

Orthonormalize given states (columns) using Gram-Schmidt.
"""
function _orthonormalize_GramSchmidt!(V::AbstractMatrix{<:Number})
    tol = 1000 * eps()
    for i in axes(V, 2)
        v = view(V, :, i) # state
        if norm(v)^2 < tol
            v .= 0 # set state to zero
            continue
        end
        for j in 1:(i - 1)
            # orthogonalize against previous state
            vj = view(V, :, j)
            a = vj ⋅ v
            axpy!(-a, vj, v) # v -= a*vj
        end
        if norm(v)^2 < tol
            v .= 0 # set state to zero
            continue
        end
        normalize!(v)
    end
    return V
end

"""
    _orthogonalize_states!(
        M1::AbstractMatrix, Q_new::AbstractMatrix, Q_old::AbstractMatrix
    )

Orthogonalize `Q_new` against `Q_old`.

Overwrites `M1`.
"""
function _orthogonalize_states!(
        M1::AbstractMatrix, Q_new::AbstractMatrix, Q_old::AbstractMatrix
    )
    mul!(M1, Q_old', Q_new)
    mul!(Q_new, Q_old, M1, -1, 1) # Q_new -= Q_old^† Q_old Q_new
    return Q_new
end
