"""
    PolesSumBlock{A <: Real, B <: Number} <: AbstractPolesSum{A, B}

Representation of block of poles on the real axis with locations ``a_i`` of type `A`
and weights ``W_i`` of type `Matrix{B}`.

```math
P(z) = ∑_i \\frac{W_i}{z-a_i}.
```

For ``W_i`` to be physical,
it needs to be Hermitian (``W_i^† = W_i``)
and positive semidefinite ``W_i ⪰ 0``.
The latter is not enforced at construction to reduce computational cost.

For a scalar variant see [`PolesSum`](@ref).
"""
struct PolesSumBlock{A <: Real, B <: Number} <: AbstractPolesSum{A, B}
    locations::Vector{A}
    weights::Vector{Matrix{B}}

    function PolesSumBlock{A, B}(locations, weights) where {A, B}
        length(locations) == length(weights) || throw(DimensionMismatch("length mismatch"))
        all(ishermitian, weights)::Bool || throw(ArgumentError("weights violate Hermiticity"))
        allequal(size, weights)::Bool ||
            throw(DimensionMismatch("weights do not have matching size"))
        _issorted_and_unique(locations)
        return new{A, B}(locations, weights)
    end
end

"""
    PolesSumBlock(locs::AbstractVector{A}, wgts::Vector{<:AbstractMatrix{B}}) where {A, B}

Create a new instance of [`PolesSumBlock`](@ref) by supplying locations `locs`
and weights `wgts`.

```jldoctest
julia> locs = 0:2;

julia> wgts = [[1 0; 0 1], [2 1; 1 2], [2 -1; -1 2]];

julia> P = PolesSumBlock(locs, wgts)
PolesSumBlock{Int64, Int64} with 3 poles of size 2×2

julia> locations(P) == locs
true

julia> weights(P) == wgts
true
```
"""
function PolesSumBlock(locs::AbstractVector{A}, wgts::Vector{<:AbstractMatrix{B}}) where {A, B}
    # Check length for permutation below.
    length(locs) == length(wgts) || throw(DimensionMismatch("length mismatch"))

    # Do no mutate user input.
    locs = collect(A, locs)
    wgts = [copy(w) for w in wgts]
    for w in wgts
        isapprox(w, w') || throw(ArgumentError("weight violates Hermiticity"))
        ishermitian(w) || hermitianpart!(w)
    end

    # sort
    p = sortperm(locs)
    locs = locs[p]
    wgts = wgts[p]

    # Merge degenerate locations.
    loc_out = similar(locs, 0)
    wgt_out = similar(wgts, 0)
    i = 1
    while i <= length(locs)
        l = locs[i]
        w = copy(wgts[i])
        i += 1
        while i <= length(locs) && locs[i] == l
            w .+= wgts[i]
            i += 1
        end
        push!(loc_out, l)
        push!(wgt_out, w)
    end
    return PolesSumBlock{A, B}(loc_out, wgt_out)
end

"""
    PolesSumBlock(
        locs::AbstractVector{A},
        amps::AbstractMatrix{B},
        tol::Real = 0
    ) where {A, B}

Create a new instance of [`PolesSumBlock`](@ref) by supplying locations `locs`
and amplitudes `amps`.

The ``i``-th column of `amps` is interpreted as the vector ``\\vec{b}_i``
and the weight as the outer product ``W_i = \\vec{b}_i \\vec{b}^†_i``.

The value `tol` determines to tolerance for location degeneracy.
If `abs(locs[i] - locs[j]) <= tol`, they are deemed to be degenerate.

```jldoctest
julia> locs = 0:1;

julia> amps = [1+2im 3im; 4 5+6im];

julia> P = PolesSumBlock(locs, amps)
PolesSumBlock{Int64, Complex{Int64}} with 2 poles of size 2×2

julia> locations(P) == locs
true

julia> weights(P) == [[5 4+8im; 4-8im 16], [9 18+15im; 18-15im 61]]
true
```
"""
function PolesSumBlock(
        locs::AbstractVector{A},
        amps::AbstractMatrix{B},
        tol::Real = 0
    ) where {A, B}
    # check input
    Base.require_one_based_indexing(locs, amps)
    tol < 0 && throw(ArgumentError("negative tol"))
    length(locs) == size(amps, 2) || throw(DimensionMismatch("locs and amps size mismatch"))
    isempty(locs) && throw(ArgumentError("no poles specified"))

    n = length(locs)

    # sort by location
    p = sortperm(locs)
    locs = locs[p]
    amps = @view amps[:, p]

    loc_new = A[]
    wgt_new = Matrix{B}[]

    i = 1
    while i <= n
        l = locs[i]
        b = @view amps[:, i]
        w = b * b'
        i += 1
        while i <= n && locs[i] - l <= tol
            # merge degenerate location
            b = @view amps[:, i]
            w .+= b * b'
            i += 1
        end
        if abs(l) <= tol
            # enforce pole at zero
            l = zero(A)
        end
        push!(loc_new, l)
        push!(wgt_new, w)
    end

    return PolesSumBlock{A, B}(loc_new, wgt_new)
end

PolesSumBlock{A, B}(P::PolesSumBlock) where {A, B} = convert(PolesSumBlock{A, B}, P)

function amplitude(P::PolesSumBlock, i::Integer, tol_amp::Real = 0; thin::Bool = false)
    tol_amp >= 0 || throw(DomainError(tol_amp, "negative amplitude"))

    w = weight(P, i)
    F = eigen(Hermitian(w))
    map!(i -> i > tol_amp^2 ? sqrt(i) : zero(i), F.values) # set small amplitudes to zero
    if !thin
        # General Julia code does not know about semipositive eigenvalues
        # and gives eltype as union of Float64 and ComplexF64.
        # Therefore, decompose by hand and apply square root in-place.
        result = F.vectors * Diagonal(F.values) * F.vectors'
        hermitianpart!(result)
    else
        n = size(w, 2)
        r = sum(>(tol_amp), F.values) # rank
        result = Matrix{eltype(F)}(undef, n, r)
        j = 1
        for i in 1:n
            F.values[i] > tol_amp || continue
            @views result[:, j] .= F.vectors[:, i] .* F.values[i]
            j += 1
        end
    end
    return result
end

function arrowhead_matrix(P::PolesSumBlock, args...; kwargs...)
    amps = amplitudes(P, args...; kwargs...)
    T = promote_type(eltype(P), eltype(eltype((amps))))
    n_b = size(P, 1)
    dim = n_b + sum(amp -> size(amp, 2), amps)
    result = zeros(T, dim, dim)::Matrix{T}

    idx = n_b
    @inbounds for i in eachindex(P)
        amp = amps[i]
        r = size(amp, 2) # rank of amplitude matrix
        result[1:n_b, (idx + 1):(idx + r)] = amp
        result[(idx + 1):(idx + r), 1:n_b] = amp'
        for j in (idx + 1):(idx + r)
            result[j, j] = location(P, i)
        end
        idx += r
    end

    return result
end

function evaluate(P::PolesSumBlock, z::Number)
    result = zeros(complex(float(eltype(P))), size(P))
    @inbounds for (ϵ, w) in P
        result .+= w .* inv(z - ϵ)
    end
    return result
end

function evaluate_gaussian(P::PolesSumBlock, ω::Real, σ::Real)
    result = zeros(complex(float(eltype(P))), size(P))
    @inbounds for (ϵ, w) in P
        re = sqrt(2) / σ * dawson((ω - ϵ) / (sqrt(2) * σ))
        im_part = pdf(Normal(ϵ, σ), ω)
        z = re - im * π * im_part
        result .+= w .* z
    end
    return result
end

function filling(P::PolesSumBlock{<:Any, B}, μ::Real = 0) where {B}
    result = zeros(B <: Real ? Float64 : ComplexF64, size(P)) # half weight changes Int → Float

    for (loc, wgt) in P
        if loc < μ
            result .+= wgt
        elseif loc == μ
            result .+= 0.5 .* wgt
        else
            break
        end
    end

    return Hermitian(result)
end

function merge_degenerate_poles!(P::PolesSumBlock, tol::Real = 0)
    # check input
    tol >= 0 || throw(ArgumentError("tol must not be negative"))
    # get information from P
    locs = locations(P)
    wgts = weights(P)

    # pole(s) at [-tol, tol]
    idx_zeros = findall(i -> abs(i) <= tol, locs)
    if !isempty(idx_zeros)
        i0 = popfirst!(idx_zeros)
        locs[i0] = 0
        for i in reverse!(idx_zeros)
            wgts[i0] .+= popat!(wgts, i)
            deleteat!(locs, i)
        end
    end

    # pole(s) at tol → ∞
    i = findfirst(>(0), locs)
    isnothing(i) && (i = lastindex(locs)) # enforce `i` to be a number
    while i < lastindex(locs)
        if locs[i + 1] - locs[i] <= tol
            # merge
            wgts[i] .+= popat!(wgts, i + 1)
            deleteat!(locs, i + 1) # keep location closer to zero
        else
            # increment index
            i += 1
        end
    end

    # pole(s) at -tol → -∞
    i = findlast(<(0), locs)
    isnothing(i) && (i = firstindex(locs)) # enforce `i` to be a number
    while i > firstindex(locs)
        if locs[i] - locs[i - 1] <= tol
            # merge
            wgts[i - 1] .+= popat!(wgts, i)
            deleteat!(locs, i - 1) # keep location closer to zero
            i -= 1
        else
            # decrement index
            i -= 1
        end
    end
    return P
end

function merge_small_weight!(P::PolesSumBlock, tol::Real)
    # check input
    tol >= 0 || throw(ArgumentError("negative tol is invalid"))
    # loop over all poles
    i = 1
    while i <= length(P)
        loc = location(P, i)
        wgt = weight(P, i)
        if norm(wgt, Inf) > tol
            # enough weight, go to next
            i += 1
            continue
        end
        if i == 1
            # add weight to next pole
            wgt_next = weight(P, i + 1)
            wgt_next .+= wgt
            deleteat!(locations(P), 1)
            deleteat!(weights(P), 1)
        elseif i == length(P)
            # add weight to previous pole
            wgt_prev = weight(P, i - 1)
            wgt_prev .+= wgt
            pop!(locations(P))
            pop!(weights(P))
        else
            # split weight such that zeroth and first moment is conserved
            loc_prev = location(P, i - 1)
            loc_next = location(P, i + 1)
            wgt_prev = weight(P, i - 1)
            wgt_next = weight(P, i + 1)
            α = (loc_next - loc) / (loc_next - loc_prev)
            wgt_prev .+= α * wgt
            wgt_next .+= (1 - α) * wgt
            deleteat!(locations(P), i)
            deleteat!(weights(P), i)
        end
    end
    return P
end

function moment(P::PolesSumBlock, n::Int = 0)
    return sum(loc^n * w for (loc, w) in P)
end

function Base.:+(A::PolesSumBlock{LA, WA}, B::PolesSumBlock{LB, WB}) where {LA, WA, LB, WB}
    L = promote_type(LA, LB)
    W = promote_type(WA, WB)
    na, nb = length(A), length(B)
    locs = Vector{L}(undef, na + nb)
    wgts = Vector{Matrix{W}}(undef, na + nb)
    ia = ib = 1
    k = 1
    @inbounds while ia <= na || ib <= nb
        if ia > na
            locs[k] = location(B, ib)
            wgts[k] = Matrix{W}(weight(B, ib))
            ib += 1
            k += 1
        elseif ib > nb
            locs[k] = location(A, ia)
            wgts[k] = Matrix{W}(weight(A, ia))
            ia += 1
            k += 1
        elseif location(A, ia) < location(B, ib)
            locs[k] = location(A, ia)
            wgts[k] = Matrix{W}(weight(A, ia))
            ia += 1
            k += 1
        elseif location(B, ib) < location(A, ia)
            locs[k] = location(B, ib)
            wgts[k] = Matrix{W}(weight(B, ib))
            ib += 1
            k += 1
        else
            # merge degenerate locations
            locs[k] = location(A, ia)
            w = Matrix{W}(weight(A, ia))
            w .+= convert(Matrix{W}, weight(B, ib))
            wgts[k] = w
            ia += 1
            ib += 1
            k += 1
        end
    end
    resize!(locs, k - 1)
    resize!(wgts, k - 1)
    return PolesSumBlock{L, W}(locs, wgts)
end

function Base.convert(::Type{PolesSumBlock{M, N}}, P::PolesSumBlock{A, B}) where {M, N, A, B}
    locs = convert(Vector{M}, locations(P))
    wgts = convert.(Matrix{N}, weights(P))
    return PolesSumBlock{M, N}(locs, wgts)
end
Base.convert(::Type{PolesSumBlock{A, B}}, P::PolesSumBlock{A, B}) where {A, B} = P

function Base.copy(P::PolesSumBlock{A, B}) where {A, B}
    return PolesSumBlock{A, B}(copy(locations(P)), map(copy, weights(P)))
end

function Base.show(io::IO, P::PolesSumBlock)
    _show_poles(io, P)
    isempty(P) || print(io, " of size ", size(P, 1), "×", size(P, 2))
    return nothing
end

Base.size(P::PolesSumBlock) = size(first(weights(P)))
Base.size(P::PolesSumBlock, i) = size(first(weights(P)), i)

function Base.transpose(P::PolesSumBlock{A, B}) where {A, B}
    wgts = [Matrix(transpose(w)) for w in weights(P)]
    return PolesSumBlock{A, B}(copy(locations(P)), wgts)
end

function LinearAlgebra.tr(P::PolesSumBlock{<:Any, B}) where {B}
    locs = copy(locations(P))
    wgts = similar(locs, real(B))

    @inbounds for i in eachindex(P)
        wgts[i] = tr(weight(P, i))
    end
    return PolesSum(locs, wgts)
end
