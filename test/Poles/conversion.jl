using RAS_DMFT
using LinearAlgebra
using Test

@testset "conversion" begin
    @testset "Anderson matrix" begin
        # scalar
        P = PolesSum([-1.0, 0.0, 1.0], [0.4, 0.2, 0.4])
        b0, HA = anderson_matrix(P)
        @test b0 ≈ 1 atol = 10 * eps()
        @test norm(
            HA - [
                0 -sqrt(0.4) -sqrt(0.4);
                -sqrt(0.4) -sqrt(0.2) 0.0;
                -sqrt(0.4) 0.0 sqrt(0.2)
            ]
        ) < 10 * eps()
        # block
        P = PolesSumBlock(
            [-1.0, 0, 1],
            [
                [1 0.1; 0.1 0.5],
                [2.0 0; 0 2],
                [1 0.1; 0.1 0.5 ],
            ]
        )
        B0, HA = anderson_matrix(P)
        @test norm(B0 * B0 - sum(weights(P))) < 50 * eps()
        F = eigen(HA)
        B = B0 * F.vectors[1:2, :]
        P_new = PolesSumBlock(copy(F.values), B)
        merge_degenerate_poles!(P_new, 50 * eps())
        @test norm(locations(P_new) - locations(P)) < 20 * eps()
        @test all(i -> i < 20 * eps(), norm.(weights(P_new) .- weights(P))) # norm of each weight difference
    end # Anderson matrix

    @testset "scalar" begin
        P = PolesContinuedFraction(5:10, 0.1:0.1:0.5, 0.25)
        PS = PolesSum(P)
        PCF = PolesContinuedFraction(PS)
        @test norm(locations(P) - locations(PCF)) < 50 * eps()
        @test norm(amplitudes(P) - amplitudes(PCF)) < 50 * eps()
        @test RAS_DMFT.scale(P) ≈ RAS_DMFT.scale(PCF) atol = 10 * eps()
        # early stopping
        PS = PolesSum([0, 1], [1, 0]) # second location has no weight
        PCF = PolesContinuedFraction(PS)
        @test length(PCF) == 1
        @test locations(PCF) == [0.0]
        @test amplitudes(PCF) == Float64[]
        @test RAS_DMFT.scale(PCF) === 1.0
    end # scalar

    @testset "block" begin
        P = PolesSumBlock([1, 2], [[1 0; 0 1], [0 0; 0 1]])
        PCF = PolesContinuedFractionBlock(P)
        @test RAS_DMFT.scale(PCF) == Diagonal([1, sqrt(2)])
        @test norm(Array(PCF) - [1 0 0 0; 0 1.5 0 0.5; 0 0 0 0; 0 0.5 0 1.5]) < 10 * eps()
        PS = PolesSumBlock(PCF)
        merge_degenerate_poles!(PS, 5 * eps())
        @test norm(locations(PS) - [1, 2]) < 10 * eps()
        @test norm(weights(PS)[1] - [1 0; 0 1]) < 10 * eps()
        @test norm(weights(PS)[2] - [0 0; 0 1]) < 10 * eps()
    end # block

    @testset "block → scalar" begin
        P = PolesSumBlock([0, 1], [[2 3im; -3im 4], [5 -6im; 6im 7]])
        foo = PolesSum(P, 1, 1)
        @test locations(foo) !== locations(P)
        @test locations(foo) == [0, 1]
        @test weights(foo) == [2, 5]
        foo = PolesSum(P, 1, 2)
        @test locations(foo) == [0, 1]
        @test weights(foo) == [3im, -6im]
        foo = PolesSum(P, 2, 1)
        @test locations(foo) == [0, 1]
        @test weights(foo) == [-3im, 6im]
        foo = PolesSum(P, 2, 2)
        @test locations(foo) == [0, 1]
        @test weights(foo) == [4, 7]
    end # block → scalar
end # conversion
