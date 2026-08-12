using RAS_DMFT
using LinearAlgebra
using Test

@testset "PolesContinuedFraction" begin
    @testset "constructor" begin
        locs = [1, 2, 3]
        amps = [4, 5]
        scl = 6

        # inner constructor
        P = PolesContinuedFraction{Int, Int}(locs, amps, scl)
        @test P.locations === locs
        @test P.amplitudes === amps
        @test P.scale === scl
        # wrong input
        @test_throws ArgumentError PolesContinuedFraction{Int, Int}(
            rand(Int, 5), rand(Int, 3), 1
        )
        @test_throws ArgumentError PolesContinuedFraction{Int, Int}(
            rand(Int, 5), rand(Int, 6), 1
        )

        # outer constructor
        P = PolesContinuedFraction(locs, amps, scl)
        @test P.locations !== locs
        @test P.locations == locs
        @test P.amplitudes !== amps
        @test P.amplitudes == amps
        @test P.scale == scl
        # default scale
        P = PolesContinuedFraction(locs, amps)
        @test P.scale === one(Int)

        # conversion of type
        P = PolesContinuedFraction(locs, amps, scl)
        P_new = PolesContinuedFraction{UInt, Float64}(P)
        @test P_new isa PolesContinuedFraction{UInt, Float64}
        @test P_new.locations == locs
        @test P_new.amplitudes == amps
        @test P_new.scale == scl
    end # constructor

    @testset "custom functions" begin
        @testset "amplitude" begin
            locs = 0:5
            amps = 6:10
            P = PolesContinuedFraction(locs, amps)
            @test amplitude(P, 1) === amps[1]
            @test amplitude(P, 5) === amps[5]
        end # amplitude

        @testset "amplitudes" begin
            locs = 0:5
            amps = 6:10
            P = PolesContinuedFraction(locs, amps)
            @test amplitudes(P) == amps
        end # amplitudes

        @testset "evaluate" begin
            locs = [-1.0, 0.0, 2.0]
            amps = [0.6, 0.4]
            scl = 1.1
            P = PolesContinuedFraction(locs, amps, scl)
            # upper/lower complex plane
            z = 0.5 + 1.0im
            @test evaluate(P, z) ≈ 0.47743341980552373 - 0.44522565484288634im atol =
                10 * eps()
            @test evaluate(P, conj(z)) == conj(evaluate(P, z))
            # Matsubara frequency
            z = 2im
            @test evaluate(P, z) ≈ 0.21044536856932053 - 0.45960359514387283im atol =
                10 * eps()
            @test evaluate(P, conj(z)) == conj(evaluate(P, z))
            # grid
            zs = [0.1 + 0.5im, 0.3 + 0.5im]
            @test evaluate(P, zs) == [evaluate(P, zs[1]), evaluate(P, zs[2])]
        end # evaluate

        @testset "evaluate_lorentzian" begin
            locs = [-1.0, 0.0, 2.0]
            amps = [0.6, 0.4]
            scl = 1.1
            # single point
            P = PolesContinuedFraction(locs, amps, scl)
            @test evaluate_lorentzian(P, 0.5, 1) ≈
                0.47743341980552373 - 0.44522565484288634im atol = 10 * eps()
            # grid
            ω = [0.1, 0.3]
            @test evaluate_lorentzian(P, ω, 0.5) ==
                [evaluate_lorentzian(P, ω[1], 0.5), evaluate_lorentzian(P, ω[2], 0.5)]
        end # evaluate_lorentzian

        @testset "locations" begin
            locs = 0:5
            amps = 6:10
            P = PolesContinuedFraction(locs, amps)
            @test locations(P) == locs
        end # locations

        @testset "scale" begin
            locs = 0:5
            amps = 6:10
            scl = 5
            P = PolesContinuedFraction(locs, amps, scl)
            @test scale(P) === scl
        end # scale

        @testset "tridiagonal_matrix" begin
            locs = 1:3
            amps = 4:5
            scl = 2
            P = PolesContinuedFraction(locs, amps, scl)
            @test tridiagonal_matrix(P) == [1 4 0; 4 2 5; 0 5 3]
            P = PolesContinuedFraction(locs, amps)
            @test tridiagonal_matrix(P) == [1 4 0; 4 2 5; 0 5 3]
        end # tridiagonal_matrix

        @testset "weight" begin
            P = PolesContinuedFraction([-1.0, 0.0, 0.5], [0.25, 1.5])
            @test_throws BoundsError weight(P, 0)
            @test weight(P, 1) == 0.0625
            @test weight(P, 2) == 2.25
            @test_throws BoundsError weight(P, 3)
        end # weight

        @testset "weights" begin
            P = PolesContinuedFraction([-1.0, 0.0, 0.5], [0.25, 1.5])
            @test weights(P) == [0.0625, 2.25]
        end # weights
    end # custom functions

    @testset "Base" begin
        @testset "eltype" begin
            locs = Float64[0.1, 0.2]
            amps = Int[1]
            P = PolesContinuedFraction(locs, amps)
            @test eltype(P) === Float64
        end # eltype

        @testset "length" begin
            locs = [0.1, 0.2]
            amps = [1]
            P = PolesContinuedFraction(locs, amps)
            @test length(P) == 2
        end # length
    end # Base

    @testset "LinearAlgebra" begin
        @testset "SymTridiagonal" begin
            locs = 1:3
            amps = 4:5
            scl = 2
            P = PolesContinuedFraction(locs, amps, scl)
            @test SymTridiagonal(P) == SymTridiagonal(locs, amps)
            P = PolesContinuedFraction(locs, amps)
            @test SymTridiagonal(P) == SymTridiagonal(locs, amps)
        end # SymTridiagonal
    end # LinearAlgebra
end # PolesContinuedFraction
