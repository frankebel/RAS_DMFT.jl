"""
    PolesSumBlock{A<:Real,B<:Number} <: AbstractPolesSum

Representation of block of poles on the real axis with locations ``a_i`` of type `A`
and weights ``W_i`` of type `Hermitian{B, Matrix{B}}` (``W_i^† = W_i``).

```math
P(ω) = ∑_i \\frac{W_i}{ω-a_i}.
```

For ``W_i`` to be physical,
it needs to be positive semidefinite ``W_i ⪰ 0``.

For a scalar variant see [`PolesSum`](@ref).
"""
struct PolesSumBlock{A <: Real, B <: Number} <: AbstractPolesSum
    locations::Vector{A}
    weights::Vector{Hermitian{B, Matrix{B}}}

    function PolesSumBlock{A, B}(locations, weights) where {A, B}
        length(locations) == length(weights) || throw(DimensionMismatch("length mismatch"))
        allequal(size, weights) ||
            throw(DimensionMismatch("weights do not have matching size"))
        result = new{A, B}(locations, weights)
        sort!(result)
        merge_degenerate_poles!(result)
        return result
    end
end

"""
    PolesSumBlock(loc::AbstractVector, wgt::Vector{<:AbstractMatrix}) where {A,B}

Create a new instance of [`PolesSumBlock`](@ref) by supplying locations `loc`
and weights `wgt`.

!!! info
    The constructor does not check if every matrix inside `wgt` is Hermitian.

```jldoctest
julia> loc = 0:2;

julia> wgt = [[1 0; 0 1], [2 1; 1 2], [2 -1; -1 2]];

julia> P = PolesSumBlock(loc, wgt)
PolesSumBlock{Int64, Int64} with 3 poles of size 2×2

julia> locations(P) == loc
true

julia> weights(P) == wgt
true
```
"""
function PolesSumBlock(loc::AbstractVector, wgt::Vector{<:AbstractMatrix})
    return PolesSumBlock(loc, map(Hermitian, wgt))
end
PolesSumBlock(loc::AbstractVector{A}, wgt::Vector{Hermitian{B, Matrix{B}}}) where {A, B} =
    PolesSumBlock{A, B}(loc, wgt)

"""
    PolesSumBlock(loc::AbstractVector, amp::AbstractMatrix{B}) where {B}

Create a new instance of [`PolesSumBlock`](@ref) by supplying locations `loc`
and amplitudes `amp`.

The ``i``-th column of `amp` is interpreted as the vector ``\\vec{b}_i``
and the weight as the outer product ``W_i = \\vec{b}_i \\vec{b}^†_i``.

```jldoctest
julia> loc = 0:1;

julia> amp = [1+2im 3im; 4 5+6im];

julia> P = PolesSumBlock(loc, amp)
PolesSumBlock{Int64, Complex{Int64}} with 2 poles of size 2×2

julia> locations(P) == loc
true

julia> weights(P) == [[5 4+8im; 4-8im 16], [9 18+15im; 18-15im 61]]
true
```
"""
function PolesSumBlock(loc::AbstractVector, amp::AbstractMatrix{B}) where {B}
    wgt = Vector{Hermitian{B, Matrix{B}}}(undef, size(amp, 2))
    for i in axes(amp, 2)
        foo = view(amp, :, i)
        wgt[i] = Hermitian(foo * foo')
    end
    return PolesSumBlock(loc, wgt)
end

PolesSumBlock{A, B}(P::PolesSumBlock) where {A, B} = convert(PolesSumBlock{A, B}, P)

# General Julia code does not know about semipositive eigenvalues
# and gives eltype as union of Float64 and ComplexF64.
# Therefore, decompose by hand and apply square root in-place.
function amplitude(P::PolesSumBlock, i::Integer)
    T = eltype(P)
    T <: Union{Float64, ComplexF64} ||
        throw(ArgumentError("amplitude(P, i) only impletmented for eltype(P) <: Union{Float64, ComplexF64}"))
    F = eigen(weight(P, i))
    tol = maximum(F.values) * sqrt(eps())
    map!(i -> i > tol ? sqrt(i) : zero(i), F.values) # set small eigenvalues to zero
    result = F.vectors * Diagonal(F.values) * F.vectors'
    return Hermitian(result)
end

function evaluate_gaussian(P::PolesSumBlock, ω::Real, σ::Real)
    d = size(P, 1)
    # dot multiplication loses Hermiticity
    real = zeros(Float64, d, d)
    imag = zero(real)
    for i in eachindex(P)
        w = weight(P, i)
        real .+= w .* sqrt(2) ./ (π * σ) .* dawson((ω - location(P, i)) / (sqrt(2) * σ))
        imag .+= w .* pdf(Normal(location(P, i), σ), ω)
    end
    result = real - im * imag
    result .*= π # not spectral function
    return result
end

function evaluate_lorentzian(P::PolesSumBlock, ω::Real, δ::Real)
    d = size(P, 1)
    result = zeros(ComplexF64, d, d)
    for i in eachindex(P)
        w = weight(P, i)
        result .+= w ./ (ω + im * δ - location(P, i))
    end
    return result
end

function merge_degenerate_poles!(P::PolesSumBlock, tol::Real = 0)
    # check input
    tol >= 0 || throw(ArgumentError("tol must not be negative"))
    # get information from P
    loc = locations(P)
    wgt = weights(P)
    # When adding matrices, Hermitian can't add in-place.
    # Thus, access H.data directly with parent(H).

    # pole(s) at [-tol, tol]
    idx_zeros = findall(i -> abs(i) <= tol, loc)
    if !isempty(idx_zeros)
        i0 = popfirst!(idx_zeros)
        loc[i0] = 0
        for i in reverse!(idx_zeros)
            parent(wgt[i0]) .+= popat!(wgt, i)
            deleteat!(loc, i)
        end
    end

    # pole(s) at tol → ∞
    i = findfirst(>(0), loc)
    isnothing(i) && (i = lastindex(loc)) # enforce `i` to be a number
    while i < lastindex(loc)
        if loc[i + 1] - loc[i] <= tol
            # merge
            parent(wgt[i]) .+= popat!(wgt, i + 1)
            deleteat!(loc, i + 1) # keep location closer to zero
        else
            # increment index
            i += 1
        end
    end

    # pole(s) at -tol → -∞
    i = findlast(<(0), loc)
    isnothing(i) && (i = firstindex(loc)) # enforce `i` to be a number
    while i > firstindex(loc)
        if loc[i] - loc[i - 1] <= tol
            # merge
            parent(wgt[i - 1]) .+= popat!(wgt, i)
            deleteat!(loc, i - 1) # keep location closer to zero
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
            parent(wgt_next) .+= wgt
            deleteat!(locations(P), 1)
            deleteat!(weights(P), 1)
        elseif i == length(P)
            # add weight to previous pole
            wgt_prev = weight(P, i - 1)
            parent(wgt_prev) .+= wgt
            pop!(locations(P))
            pop!(weights(P))
        else
            # split weight such that zeroth and first moment is conserved
            loc_prev = location(P, i - 1)
            loc_next = location(P, i + 1)
            wgt_prev = weight(P, i - 1)
            wgt_next = weight(P, i + 1)
            α = (loc_next - loc) / (loc_next - loc_prev)
            parent(wgt_prev) .+= α * wgt
            parent(wgt_next) .+= (1 - α) * wgt
            deleteat!(locations(P), i)
            deleteat!(weights(P), i)
        end
    end
    return P
end

function moment(P::PolesSumBlock, n::Int = 0)
    return sum(i -> i[1]^n * i[2], zip(locations(P), weights(P)))
end

function Core.Array(P::PolesSumBlock)
    T = eltype(P) <: Real ? Float64 : ComplexF64
    N = length(P)
    n = size(P, 1)
    result = zeros(T, (N + 1) * n, (N + 1) * n)
    A = locations(P)
    result[1:n, 1:n] = zeros(T, n, n)
    for i in 2:(N + 1)
        i1 = (i - 1) * n + 1
        i2 = i * n
        B = amplitude(P, i - 1)
        result[1:n, i1:i2] = B # first block row
        result[i1:i2, 1:n] = B # first block column
        map(j -> result[j, j] = A[i - 1], i1:i2) # main diagonal
    end
    return result
end

function Base.:+(A::PolesSumBlock{<:Any, TA}, B::PolesSumBlock{<:Any, TB}) where {TA, TB}
    loc = [locations(A); locations(B)]
    # copy weights of `A`, `B` to `wgt`
    T = promote_type(TA, TB)
    wgt = Vector{Hermitian{T, Matrix{T}}}(undef, length(A) + length(B))
    for i in eachindex(wgt)
        if i <= length(A)
            wgt[i] = copy(weight(A, i))
        else
            wgt[i] = copy(weight(B, i - length(A)))
        end
    end
    return PolesSumBlock(loc, wgt)
end

function Base.convert(::Type{PolesSumBlock{M, N}}, P::PolesSumBlock{A, B}) where {M, N, A, B}
    loc = convert(Vector{M}, locations(P))
    wgt = convert.(Hermitian{N, Matrix{N}}, weights(P))
    return PolesSumBlock(loc, wgt)
end
Base.convert(::Type{PolesSumBlock{A, B}}, P::PolesSumBlock{A, B}) where {A, B} = P

function Base.copy(P::PolesSumBlock)
    return PolesSumBlock(copy(locations(P)), map(copy, weights(P)))
end

Base.eltype(::Type{<:PolesSumBlock{A, B}}) where {A, B} = promote_type(A, B)

function Base.show(io::IO, P::PolesSumBlock)
    print(
        io, summary(P), " with ", length(P), " poles"
    )
    isempty(P) || print(io, " of size ", size(P, 1), "×", size(P, 2))
    return Nothing
end

Base.size(P::PolesSumBlock) = size(first(weights(P)))
Base.size(P::PolesSumBlock, i) = size(first(weights(P)), i)

function Base.transpose(P::PolesSumBlock)
    return PolesSumBlock(copy(locations(P)), map(conj, weights(P)))
end
