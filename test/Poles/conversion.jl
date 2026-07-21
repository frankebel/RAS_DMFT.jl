using RAS_DMFT
using LinearAlgebra
using Test

@testset "conversion" begin
    @testset "Anderson matrix" begin
        # scalar
        P = PolesSum([-1.0, 0.0, 1.0], [0.4, 0.2, 0.4])
        b_0, H_A = anderson_matrix(P)
        @test b_0 ≈ 1 atol = 10 * eps()
        @test norm(
            H_A - [
                0 -sqrt(0.4) sqrt(0.4);
                -sqrt(0.4) -sqrt(0.2) 0.0;
                sqrt(0.4) 0.0 sqrt(0.2)
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
        B_0, H_A = anderson_matrix(P)
        @test norm(B_0 * B_0 - sum(weights(P))) < 50 * eps()
        F = eigen(H_A)
        B = B_0 * F.vectors[1:2, :]
        P_new = PolesSumBlock(copy(F.values), B, 50 * eps())
        @test norm(locations(P_new) - locations(P)) < 20 * eps()
        @test all(<(50 * eps()), norm.(weights(P_new) .- weights(P))) # norm of each weight difference
    end # Anderson matrix

    @testset "arrowhead matrix" begin
        # scalar
        loc = 1:5
        amp = 6:10
        P = PolesSum(loc, abs2.(amp))
        ma = arrowhead_matrix(P)
        @inferred arrowhead_matrix(P)
        @test ma isa Matrix{Float64}
        @test ma == [
            0 6 7 8 9 10
            6 1 0 0 0 0
            7 0 2 0 0 0
            8 0 0 3 0 0
            9 0 0 0 4 0
            10 0 0 0 0 5
        ]
        # pole with zero weight
        loc = 1:5
        amp = [6, 7, 0, 9, 0]
        P = PolesSum(loc, abs2.(amp))
        ma = arrowhead_matrix(P)
        @test ma == [
            0 6 7 0 9 0
            6 1 0 0 0 0
            7 0 2 0 0 0
            0 0 0 3 0 0
            9 0 0 0 4 0
            0 0 0 0 0 5
        ]

        # block
        locs = 1:2
        amps = [
            rand(ComplexF64, 4, 1),
            rand(ComplexF64, 4, 2),
        ]
        wgts = map(amp -> amp * amp', amps)
        P = PolesSumBlock(locs, wgts)
        ma = arrowhead_matrix(P)
        @inferred arrowhead_matrix(P)
        @test ishermitian(ma)
        @test iszero(@view ma[1:4, 1:4])
        amp1 = @view ma[1:4, 5:8]
        @test norm(wgts[1] - amp1 * amp1') < 100 * eps()
        amp2 = @view ma[1:4, 9:12]
        @test norm(wgts[2] - amp2 * amp2') < 100 * eps()
        @test view(ma, 5:12, 5:12) == Diagonal([1, 1, 1, 1, 2, 2, 2, 2])
        # thin rectangular amplitudes
        ma = arrowhead_matrix(P, 10 * sqrt(eps()); thin = true)
        @test ishermitian(ma)
        @test iszero(@view ma[1:4, 1:4])
        amp1 = @view ma[1:4, 5]
        @test norm(wgts[1] - amp1 * amp1') < 100 * eps()
        amp2 = @view ma[1:4, 6:7]
        @test norm(wgts[2] - amp2 * amp2') < 100 * eps()
        @test view(ma, 5:7, 5:7) == Diagonal([1, 2, 2])
    end # arrowhead matrix

    @testset "scalar" begin
        P = PolesContinuedFraction(5:10, 0.1:0.1:0.5, 0.25)
        PS = PolesSum(P)
        PCF = PolesContinuedFraction(PS)
        @test norm(locations(P) - locations(PCF)) < 50 * eps()
        @test norm(amplitudes(P) - amplitudes(PCF)) < 50 * eps()
        @test scale(P) ≈ scale(PCF) atol = 10 * eps()
        # early stopping
        PS = PolesSum([0, 1], [1, 0]) # second location has no weight
        PCF = PolesContinuedFraction(PS)
        @test length(PCF) == 1
        @test locations(PCF) == [0.0]
        @test amplitudes(PCF) == Float64[]
        @test scale(PCF) === 1.0
    end # scalar

    @testset "block" begin
        P = PolesSumBlock([1, 2], [[1.0 0; 0 1], [0 0; 0 1]])
        PCF = PolesContinuedFractionBlock(P)
        @test scale(PCF) == Diagonal([1, sqrt(2)])
        @test norm(tridiagonal_matrix(PCF) - [1 0 0 0; 0 1.5 0 0.5; 0 0 0 0; 0 0.5 0 1.5]) < 10 * eps()
        PS = PolesSumBlock(PCF)
        merge_degenerate_poles!(PS, 5 * eps())
        @test norm(locations(PS) - [1, 2]) < 10 * eps()
        @test norm(weight(PS, 1) - [1 0; 0 1]) < 10 * eps()
        @test norm(weight(PS, 2) - [0 0; 0 1]) < 10 * eps()
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
