using RAS_DMFT
using Test

@testset "IO" begin
    @testset "Number" begin
        s = 1 + 2im
        @test write_hdf5("test.h5", s) === nothing
        @inferred read_hdf5("test.h5", Complex{Int})
        @test read_hdf5("test.h5", Complex{Int}) === 1 + 2im
        @test read_hdf5("test.h5", ComplexF64) === 1.0 + 2.0im
        @test_throws InexactError read_hdf5("test.h5", Int)
    end # Number

    @testset "Array" begin
        # Vector{Number}
        v = rand(10)
        @test write_hdf5("test.h5", v) === nothing
        @inferred read_hdf5("test.h5", Vector{Float64})
        @test read_hdf5("test.h5", Vector{Float64}) == v

        # Matrix{Number}
        m = rand(10, 10)
        @test write_hdf5("test.h5", m) === nothing
        @inferred read_hdf5("test.h5", Matrix{Float64})
        @test read_hdf5("test.h5", Matrix{Float64}) == m

        # Vector{Matrix{Number}}
        v = [rand(2, 2) for _ in 1:10]
        @test write_hdf5("test.h5", v) === nothing
        @inferred read_hdf5("test.h5", Vector{Matrix{Float64}})
        @test read_hdf5("test.h5", Vector{Matrix{Float64}}) == v
    end # Array

    @testset "PolesSum" begin
        P = PolesSum([1.0, 2], [3, 4])
        # write
        @test write_hdf5("test.h5", P) === nothing
        # read
        @inferred read_hdf5("test.h5", PolesSum{Float64, Int})
        foo = read_hdf5("test.h5", PolesSum{Float64, Int})
        @test locations(foo) == [1, 2]
        @test weights(foo) == [3.0, 4.0]
    end # PolesSum

    @testset "PolesSumBlock" begin
        P = PolesSumBlock([1.0, 2], [[1 0; 0 1], [2 1; 1 0]])
        # write
        @test write_hdf5("test.h5", P) === nothing
        # read
        @inferred read_hdf5("test.h5", PolesSumBlock{Float64, Int})
        foo = read_hdf5("test.h5", PolesSumBlock{Float64, Int})
        @test locations(foo) == [1.0, 2]
        @test weights(foo) == [[1 0; 0 1], [2 1; 1 0]]
    end # PolesSumBlock

    @testset "PolesContinuedFraction" begin
        P = PolesContinuedFraction([1.0, 2], [3], 4)
        # write
        @test write_hdf5("test.h5", P) === nothing
        # read
        @inferred read_hdf5("test.h5", PolesContinuedFraction{Float64, Int})
        foo = read_hdf5("test.h5", PolesContinuedFraction{Float64, Int})
        @test locations(foo) == [1.0, 2]
        @test amplitudes(foo) == [3]
        @test scale(foo) == 4
    end # PolesContinuedFraction

    @testset "PolesContinuedFractionBlock" begin
        P = PolesContinuedFractionBlock([[1 0; 0 1], [2 1; 1 0]], [[3 1; 1 0]], [0 0; 0 0])
        # write
        @test write_hdf5("test.h5", P) === nothing
        # read
        @inferred read_hdf5("test.h5", PolesContinuedFractionBlock{Int, Int})
        foo = read_hdf5("test.h5", PolesContinuedFractionBlock{Int, Int})
        @test locations(foo) == [[1 0; 0 1], [2 1; 1 0]]
        @test amplitudes(foo) == [[3 1; 1 0]]
        @test scale(foo) == [0 0; 0 0]
    end # PolesContinuedFractionBlock

    isfile("test.h5") && rm("test.h5")
end # IO
