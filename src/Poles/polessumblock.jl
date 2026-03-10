"""
    PolesSumBlock{A<:Real,B<:Number} <: AbstractPolesSum

Representation of block of poles on the real axis with locations ``a_i`` of type `A`
and weights ``B_i`` of type `Matrix{B}`.

```math
P(ω) = ∑_i \\frac{B_i}{ω-a_i}.
```

For a scalar variant see [`PolesSum`](@ref).
"""
struct PolesSumBlock{A <: Real, B <: Number} <: AbstractPolesSum
    locations::Vector{A}
    weights::Vector{Matrix{B}}

    function PolesSumBlock{A, B}(locations, weights) where {A, B}
        length(locations) == length(weights) || throw(DimensionMismatch("length mismatch"))
        all(ishermitian, weights) || throw(ArgumentError("weights are not hermitian"))
        allequal(size, weights) ||
            throw(DimensionMismatch("weights do not have matching size"))
        return new{A, B}(locations, weights)
    end
end

"""
    PolesSumBlock(loc::AbstractVector{A}, wgt::Vector{<:AbstractMatrix{B}}) where {A,B}

Create a new instance of [`PolesSumBlock`](@ref) by supplying locations `loc`
and weights `wgt`.

```jldoctest
julia> loc = collect(0:2);

julia> wgt = [[1 0; 0 1], [2 1; 1 2], [2 -1; -1 2]];

julia> P = PolesSumBlock(loc, wgt)
PolesSumBlock{Int64, Int64} with 3 poles of size 2×2

julia> locations(P) === loc
true

julia> weights(P) === wgt
true
```
"""
PolesSumBlock(loc::AbstractVector{A}, wgt::Vector{<:AbstractMatrix{B}}) where {A, B} =
    PolesSumBlock{A, B}(loc, wgt)

"""
    PolesSumBlock(loc::AbstractVector, amp::AbstractMatrix{B}) where {B}

Create a new instance of [`PolesSumBlock`](@ref) by supplying locations `loc`
and amplitudes `amp`.

The ``i``-th column of `amp` is interpreted as the vector ``\\vec{b}_i`` and the weight
as ``B_i = \\vec{b}_i \\vec{b}^†_i``.

```jldoctest
julia> loc = collect(0:1);

julia> amp = [1+2im 3im; 4 5+6im];

julia> P = PolesSumBlock(loc, amp)
PolesSumBlock{Int64, Complex{Int64}} with 2 poles of size 2×2

julia> locations(P) === loc
true

julia> weights(P) == [[5 4+8im; 4-8im 16], [9 18+15im; 18-15im 61]]
true
```
"""
function PolesSumBlock(loc::AbstractVector, amp::AbstractMatrix{B}) where {B}
    wgt = Vector{Matrix{B}}(undef, size(amp, 2))
    for i in axes(amp, 2)
        foo = view(amp, :, i)
        wgt[i] = foo * foo'
    end
    return PolesSumBlock(loc, wgt)
end

function PolesSumBlock{A, B}(P::PolesSumBlock) where {A, B}
    return PolesSumBlock(Vector{A}(locations(P)), map(i -> Matrix{B}(i), weights(P)))
end

# `sqrt(matrix)` would work but is very slow on first run.
# Annotating as symmetric or hermitian makes it much faster.
# But general Julia code still does not know about semipositive eigenvalues
# and gives eltype as union of Float64 and ComplexF64.
# Therefore, decompose by hand and apply square root in-place.
function amplitude(P::PolesSumBlock, i::Integer)
    T = eltype(P) <: Real ? Float64 : ComplexF64 # use double precision
    F = eigen!(Hermitian(Matrix{T}(weight(P, i))))
    tol = maximum(F.values) * sqrt(eps())
    map!(i -> i > tol ? sqrt(i) : zero(i), F.values) # set small eigenvalues to zero
    result = F.vectors * Diagonal(F.values) * F.vectors'
    hermitianpart!(result)
    return result
end

function evaluate_gaussian(P::PolesSumBlock, ω::Real, σ::Real)
    d = size(P, 1)
    real = zeros(Float64, d, d)
    imag = zero(real)
    for i in eachindex(P)
        w = weight(P, i)
        real .+= w .* sqrt(2) ./ (π * σ) .* dawson((ω - locations(P)[i]) / (sqrt(2) * σ))
        imag .+= w .* pdf(Normal(locations(P)[i], σ), ω)
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
        result .+= w ./ (ω + im * δ - locations(P)[i])
    end
    return result
end

function merge_degenerate_poles!(P::PolesSumBlock, tol::Real = 0)
    # check input
    tol >= 0 || throw(ArgumentError("tol must not be negative"))
    issorted(P) || throw(ArgumentError("P must be sorted"))
    # get information from P
    loc = locations(P)
    wgt = weights(P)
    # pole(s) at [-tol, tol]
    idx_zeros = findall(i -> abs(i) <= tol, loc)
    if !isempty(idx_zeros)
        i0 = popfirst!(idx_zeros)
        loc[i0] = 0
        for i in reverse!(idx_zeros)
            wgt[i0] .+= popat!(wgt, i)
            deleteat!(loc, i)
        end
    end
    # pole(s) at tol → ∞
    i = findfirst(>(0), loc)
    isnothing(i) && (i = lastindex(loc)) # enforce `i` to be a number
    while i < lastindex(loc)
        if loc[i + 1] - loc[i] <= tol
            # merge
            wgt[i] .+= popat!(wgt, i + 1)
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
            wgt[i - 1] .+= popat!(wgt, i)
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
    issorted(P) || throw(ArgumentError("P must be sorted"))
    # loop over all poles
    i = 1
    while i <= length(P)
        loc = locations(P)[i]
        wgt = weight(P, i)
        if norm(wgt, Inf) > tol
            # enough weight, go to next
            i += 1
            continue
        end
        if i == 1
            # add weight to next pole
            wgt_next = weight(P, i + 1)
            @. wgt_next += wgt
            deleteat!(locations(P), 1)
            deleteat!(weights(P), 1)
        elseif i == length(P)
            # add weight to previous pole
            wgt_prev = weight(P, i - 1)
            @. wgt_prev += wgt
            pop!(locations(P))
            pop!(weights(P))
        else
            # split weight such that zeroth and first moment is conserved
            loc_prev = locations(P)[i - 1]
            loc_next = locations(P)[i + 1]
            wgt_prev = weight(P, i - 1)
            wgt_next = weight(P, i + 1)
            α = (loc_next - loc) / (loc_next - loc_prev)
            @. wgt_prev += α * wgt
            @. wgt_next += (1 - α) * wgt
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
    wgt = Vector{Matrix{promote_type(TA, TB)}}(undef, length(A) + length(B))
    for i in eachindex(wgt)
        if i <= length(A)
            wgt[i] = copy(weight(A, i))
        else
            wgt[i] = copy(weight(B, i - length(A)))
        end
    end
    result = PolesSumBlock(loc, wgt)
    sort!(result)
    merge_degenerate_poles!(result, 0)
    return result
end

function Base.copy(P::PolesSumBlock)
    return PolesSumBlock(copy(locations(P)), map(i -> copy(i), weights(P)))
end

Base.eltype(::Type{<:PolesSumBlock{A, B}}) where {A, B} = promote_type(A, B)

function Base.show(io::IO, P::PolesSumBlock)
    return print(
        io, summary(P), " with ", length(P), " poles of size ", size(P, 1), "×", size(P, 2)
    )
end

Base.size(P::PolesSumBlock) = size(first(weights(P)))
Base.size(P::PolesSumBlock, i) = size(first(weights(P)), i)

function Base.transpose(P::PolesSumBlock)
    return PolesSumBlock(locations(P), map(transpose, weights(P)))
end
