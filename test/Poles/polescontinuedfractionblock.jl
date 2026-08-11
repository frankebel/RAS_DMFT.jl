using RAS_DMFT
using LinearAlgebra
using Test

@testset "PolesContinuedFractionBlock" begin
    @testset "constructor" begin
        locs = [[1 2; 2 3], [4 5; 5 6]]
        amps = [[7 8; 8 9]]
        scl = [10 11; 11 12]

        # inner constructor
        P = PolesContinuedFractionBlock{Int, Int}(locs, amps, scl)
        @test P.locations === locs
        @test P.amplitudes === amps
        @test P.scale === scl
        # wrong input
        # matrices not hermitian
        @test_throws ArgumentError PolesContinuedFractionBlock{Int, Int}(locs, locs, scl)
        @test_throws ArgumentError PolesContinuedFractionBlock{Int, Int}(
            [[1 2; 3 4], [4 5; 5 6]], amps, scl
        )
        @test_throws ArgumentError PolesContinuedFractionBlock{Int, Int}(
            locs, [[7 8; 9 10]], scl
        )
        @test_throws ArgumentError PolesContinuedFractionBlock{Int, Int}(
            locs, amps, [10 11; 12 13]
        )
        # matrices have wrong size
        @test_throws DimensionMismatch PolesContinuedFractionBlock{Int, Int}(
            [[1;;], [4 5; 5 6]], amps, scl
        )
        @test_throws DimensionMismatch PolesContinuedFractionBlock{Int, Int}(
            locs, [[1;;]], scl
        )
        @test_throws DimensionMismatch PolesContinuedFractionBlock{Int, Int}(locs, amps, [1;;])

        # outer constructor
        P = PolesContinuedFractionBlock(locs, amps, scl)
        @test P.locations === locs
        @test P.amplitudes === amps
        @test P.scale === scl
        # # default scale
        P = PolesContinuedFractionBlock(locs, amps)
        @test P.locations === locs
        @test P.amplitudes === amps
        @test P.scale == [1 0; 0 1]

        # conversion of type
        P = PolesContinuedFractionBlock(locs, amps, scl)
        P_new = PolesContinuedFractionBlock{UInt, Float64}(P)
        @test typeof(P_new) === PolesContinuedFractionBlock{UInt, Float64}
        @test P_new.locations == locs
        @test P_new.amplitudes == amps
        @test P_new.scale == scl
    end # constructor

    @testset "custom functions" begin
        @testset "amplitude" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            P = PolesContinuedFractionBlock(locs, amps)
            @test amplitude(P, 1) === amps[1]
        end # amplitude

        @testset "amplitudes" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            P = PolesContinuedFractionBlock(locs, amps)
            @test amplitudes(P) === amps
        end # amplitudes

        @testset "evaluate" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            scl = [10 11; 11 12]
            P = PolesContinuedFractionBlock(locs, amps, scl)
            # upper/lower complex plane
            z = 0.5 + 1.0im
            @test norm(
                evaluate(P, z) - [
                    9.3868717682668326 - 1.7247467522667106im 10.275713626756321 - 1.872408526521784im
                    10.275713626756325 - 1.8724085265217789im 11.250661306616291 - 2.0359262909579234im
                ],
            ) < 100 * eps()
            @test norm(evaluate(P, conj(z)) - conj.(permutedims(evaluate(P, z)))) <
                100 * eps()
            # Matsubara frequency
            z = 2im
            @test norm(
                evaluate(P, z) - [
                    9.6216882806905861 - 3.3491825125308288im 10.540934044076694 - 3.6580774126819939im
                    10.540934044076694 - 3.6580774126819948im 11.548094518259207 - 3.9976827909937138im
                ],
            ) < 100 * eps()
            @test norm(evaluate(P, conj(z)) - conj.(permutedims(evaluate(P, z)))) <
                100 * eps()
            # grid
            zs = [0.1 + 0.5im, 0.3 + 0.5im]
            @test evaluate(P, zs) == [evaluate(P, zs[1]), evaluate(P, zs[2])]
        end # evaluate

        @testset "evaluate_lorentzian" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            scl = [10 11; 11 12]
            P = PolesContinuedFractionBlock(locs, amps, scl)
            # single point
            @test norm(
                evaluate_lorentzian(P, 0.5, 1) - [
                    9.3868717682668326 - 1.7247467522667106im 10.275713626756321 - 1.872408526521784im
                    10.275713626756325 - 1.8724085265217789im 11.250661306616291 - 2.0359262909579234im
                ],
            ) < 10 * eps()
            # grid
            ω = [0.1, 0.3]
            @test evaluate_lorentzian(P, ω, 0.5) ==
                [evaluate_lorentzian(P, ω[1], 0.5), evaluate_lorentzian(P, ω[2], 0.5)]
        end # evaluate_lorentzian

        @testset "locations" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            P = PolesContinuedFractionBlock(locs, amps)
            @test locations(P) === locs
        end # locations

        @testset "scale" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            scl = [10 11; 11 12]
            P = PolesContinuedFractionBlock(locs, amps, scl)
            @test scale(P) === scl
        end # scale

        @testset "tridiagonal_matrix" begin
            locs = [[1 2; 2 3], [4 5; 5 6], [7 8; 8 9]]
            amps = [[10 11; 11 12], [13 14; 14 15]]
            scl = [16 17; 17 18]
            P = PolesContinuedFractionBlock(locs, amps, scl)
            @test tridiagonal_matrix(P) == [
                1  2  10 11 0  0
                2  3  11 12 0  0
                10 11 4  5  13 14
                11 12 5  6  14 15
                0  0  13 14 7  8
                0  0  14 15 8  9
            ]
            P = PolesContinuedFractionBlock(locs, amps)
            @test tridiagonal_matrix(P) == [
                1  2  10 11 0  0
                2  3  11 12 0  0
                10 11 4  5  13 14
                11 12 5  6  14 15
                0  0  13 14 7  8
                0  0  14 15 8  9
            ]
        end # tridiagonal_matrix

        @testset "weight" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            scl = [10 11; 11 12]
            P = PolesContinuedFractionBlock(locs, amps, scl)
            @test_throws MethodError weight(P, 1)
        end # weight

        @testset "weights" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            scl = [10 11; 11 12]
            P = PolesContinuedFractionBlock(locs, amps, scl)
            @test_throws MethodError weights(P)
        end # weights
    end # custom functions

    @testset "Base" begin
        @testset "eltype" begin
            locs = Matrix{Float64}[[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            P = PolesContinuedFractionBlock(locs, amps)
            @test eltype(P) === Float64
        end # eltype

        @testset "length" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            P = PolesContinuedFractionBlock(locs, amps)
            @test length(P) == 2
        end # length

        @testset "size" begin
            locs = [[1 2; 2 3], [4 5; 5 6]]
            amps = [[7 8; 8 9]]
            P = PolesContinuedFractionBlock(locs, amps)
            @test size(P) == (2, 2)
            @test size(P, 1) === 2
            @test size(P, 2) === 2
        end # size
    end # Base
end # PolesContinuedFractionBlock
