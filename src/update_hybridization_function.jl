# Bethe lattice trick: Δ(z) = Δ0(z + μ - Σ(z))

"""
    update_hybridization_function(
        Δ0::PolesSum{R,R}, μ::R, Σ_H::R, Σ::PolesSum{R,R}
    ) where {R<:Real}


Calculate the new hybridization function in [`PolesSum`](@ref) representation.

```math
Δ(ω) = Δ_0(ω + μ - Σ(ω))
```
"""
function update_hybridization_function(
        Δ0::PolesSum{R, R}, μ::R, Σ_H::R, Σ::PolesSum{R, R}
    ) where {R <: Real}
    Σ = remove_zero_weight(Σ)

    n = length(Σ) + 1
    n_tot = length(Δ0) * n
    loc_new = Vector{R}(undef, n_tot)
    wgt_new = Vector{R}(undef, n_tot)

    # Create tridiagonal matrix `T` using Householder transformations.
    # These do not touch the (1,1) element of the original matrix,
    # making it perfect for our use case.
    h = hessenberg!(Array(Σ)) # don't give symmetric information on purpose
    T = SymTridiagonal(diag(h.H), diag(h.H, -1)) # diagonal and first lower diagonal

    # Update the (1,1) element for each pole in Δ0 and diagonalize `T`.
    Threads.@threads for i in eachindex(Δ0)
        idx_low = 1 + n * (i - 1)
        idx_high = idx_low + n - 1
        bar = copy(T)
        bar[1, 1] = Σ_H - μ + locations(Δ0)[i]
        loc_new[idx_low:idx_high], U = eigen!(bar)
        U = h.Q * U # transform back
        wgt_new[idx_low:idx_high] = weight(Δ0, i) .* abs2.(view(U, 1, :)) # multiply new weights with original
    end

    return PolesSum(loc_new, wgt_new)
end
