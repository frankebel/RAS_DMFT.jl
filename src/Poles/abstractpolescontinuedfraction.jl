"""
    AbstractPolesContinuedFraction{A, B} <: AbstractPoles{A, B}

Supertype which represents a (block) function on the real axis as a continued fraction.
"""
abstract type AbstractPolesContinuedFraction{A, B} <: AbstractPoles{A, B} end

amplitude(P::AbstractPolesContinuedFraction, i::Integer) = amplitudes(P)[i]

amplitudes(P::AbstractPolesContinuedFraction) = P.amplitudes

"""
    scale(P::AbstractPolesContinuedFraction)

Return the scale of `P`.
"""
scale(P::AbstractPolesContinuedFraction) = P.scale

"""
    tridiagonal_matrix(P::AbstractPolesContinuedFraction)

Calculate the (block) tridiagonal representation of `P`

```math
\\begin{pmatrix}
A_1 & B_1 &     &         &         \\\\
B_1 & A_2 & B_2 &         &         \\\\
    & B_2 & ⋱   &  ⋱      &         \\\\
    &     & ⋱   &  ⋱      & B_{N-1} \\\\
    &     &     & B_{N-1} & A_N
\\end{pmatrix} .
```

The overall scale of `P` is not part of the matrix.
"""
function tridiagonal_matrix end
