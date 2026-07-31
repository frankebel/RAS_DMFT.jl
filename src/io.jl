# IO to read/write data in HDF5 format.

"""
    read_hdf5(filename::Abstractstring, T::Type)

Read data of type `T` from `filename`.

See also [`write_hdf5`](@ref).
"""
function read_hdf5 end

"""
    write_hdf5(filename::Abstractstring, content)

Write `content` to `filename` in HDF5 format.

If the file already exists, it will be overwritten.

See also [`read_hdf5`](@ref).
"""
function write_hdf5 end

# Number
function read_hdf5(filename::AbstractString, ::Type{T}) where {T <: Number}
    return h5open(filename, "r") do fid
        scalar::T = read(fid, "s")
        return scalar
    end
end

function write_hdf5(filename::AbstractString, s::Number)
    h5open(filename, "w") do fid
        fid["s"] = s
    end
    return nothing
end

# Array{Number,N}
function read_hdf5(filename::AbstractString, ::Type{<:Array{T, N}}) where {T <: Number, N}
    return h5open(filename, "r") do fid
        a::Array{T, N} = read(fid, "a")
        return a
    end
end

function write_hdf5(filename::AbstractString, content::Array{T, N}) where {T, N}
    h5open(filename, "w") do fid
        fid["a"] = content
    end
    return nothing
end

# Vector{Matrix{Number}}
function read_hdf5(
        filename::AbstractString, ::Type{<:Vector{<:Matrix{<:T}}}
    ) where {T <: Number}
    return h5open(filename, "r") do fid
        n::Int = read(fid, "n") # read out length
        result = Vector{Matrix{T}}(undef, n)
        for i in 1:n
            result[i] = read(fid, "$i")
        end
        return result
    end
end

function write_hdf5(filename::AbstractString, content::Vector{Matrix{T}}) where {T <: Number}
    n = length(content)
    h5open(filename, "w") do fid
        fid["n"] = n # store length
        for i in eachindex(content)
            fid["$i"] = content[i] # store each matrix
        end
    end
    return nothing
end

# PolesSum{A,B}
function read_hdf5(filename::AbstractString, ::Type{<:PolesSum{A, B}}) where {A, B}
    return h5open(filename, "r") do fid
        locations::Vector{A} = read(fid, "locations")
        weights::Vector{B} = read(fid, "weights")
        return PolesSum(locations, weights)
    end
end

function write_hdf5(filename::AbstractString, P::PolesSum)
    h5open(filename, "w") do fid
        fid["locations"] = locations(P)
        fid["weights"] = weights(P)
    end
    return nothing
end

# PolesSumBlock{A,B}
function read_hdf5(filename::AbstractString, ::Type{<:PolesSumBlock{A, B}}) where {A, B}
    return h5open(filename, "r") do fid
        locations::Vector{A} = read(fid, "locations")
        weights = Vector{Matrix{B}}(undef, length(locations))
        for i in eachindex(locations)
            weights[i] = read(fid, "weights/$i")
        end
        return PolesSumBlock(locations, weights)
    end
end

function write_hdf5(filename::AbstractString, P::PolesSumBlock)
    h5open(filename, "w") do fid
        fid["locations"] = locations(P)
        for i in eachindex(P)
            fid["weights/$i"] = weight(P, i)
        end
    end
    return nothing
end

# PolesContinuedFraction{A,B}
function read_hdf5(
        filename::AbstractString, ::Type{<:PolesContinuedFraction{A, B}}
    ) where {A, B}
    return h5open(filename, "r") do fid
        locations::Vector{A} = read(fid, "locations")
        amplitudes::Vector{B} = read(fid, "amplitudes")
        scale::B = read(fid, "scale")
        return PolesContinuedFraction(locations, amplitudes, scale)
    end
end

function write_hdf5(filename::AbstractString, P::PolesContinuedFraction)
    h5open(filename, "w") do fid
        fid["locations"] = locations(P)
        fid["amplitudes"] = amplitudes(P)
        fid["scale"] = scale(P)
    end
    return nothing
end

# PolesContinuedFractionBlock{A,B}
function read_hdf5(
        filename::AbstractString, ::Type{<:PolesContinuedFractionBlock{A, B}}
    ) where {A, B}
    return h5open(filename, "r") do fid
        n::Int = read(fid, "length")
        locations = Vector{Matrix{A}}(undef, n)
        amplitudes = Vector{Matrix{B}}(undef, n - 1)
        for i in 1:(n - 1)
            locations[i] = read(fid, "locations/$i")
            amplitudes[i] = read(fid, "amplitudes/$i")
        end
        locations[n] = read(fid, "locations/$n")
        scale::Matrix{B} = read(fid, "scale")
        return PolesContinuedFractionBlock(locations, amplitudes, scale)
    end
end

function write_hdf5(filename::AbstractString, P::PolesContinuedFractionBlock)
    h5open(filename, "w") do fid
        n = length(P)
        fid["length"] = n
        for i in 1:(n - 1)
            fid["locations/$i"] = locations(P)[i]
            fid["amplitudes/$i"] = amplitude(P, i)
        end
        fid["locations/$n"] = locations(P)[n]
        fid["scale"] = scale(P)
    end
    return nothing
end
