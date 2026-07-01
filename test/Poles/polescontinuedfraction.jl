using RAS_DMFT
using LinearAlgebra
using Test

@testset "PolesContinuedFraction" begin
    @testset "constructor" begin
        loc = 0:5
        amp = 6:10
        scl = 11

        # inner constructor
        P = PolesContinuedFraction{Int, Int}(loc, amp, scl)
        @test P.locations == loc
        @test P.amplitudes == amp
        @test P.scale == scl
        # wrong input
        @test_throws ArgumentError PolesContinuedFraction{Int, Int}(
            rand(Int, 5), rand(Int, 3), 1
        )
        @test_throws ArgumentError PolesContinuedFraction{Int, Int}(
            rand(Int, 5), rand(Int, 6), 1
        )

        # outer constructor
        P = PolesContinuedFraction(loc, amp, scl)
        @test P.locations == loc
        @test P.amplitudes == amp
        @test P.scale == scl
        # default scale
        P = PolesContinuedFraction(loc, amp)
        @test P.locations == loc
        @test P.amplitudes == amp
        @test P.scale == one(Int)

        # conversion of type
        P = PolesContinuedFraction(loc, amp, scl)
        P_new = PolesContinuedFraction{UInt, Float64}(P)
        @test typeof(P_new) === PolesContinuedFraction{UInt, Float64}
        @test P_new.locations == loc
        @test P_new.amplitudes == amp
        @test P_new.scale == scl
    end # constructor

    @testset "custom functions" begin
        @testset "amplitude" begin
            loc = 0:5
            amp = 6:10
            P = PolesContinuedFraction(loc, amp)
            @test amplitude(P, 1) === amp[1]
            @test amplitude(P, 5) === amp[5]
        end # amplitude

        @testset "amplitudes" begin
            loc = 0:5
            amp = 6:10
            P = PolesContinuedFraction(loc, amp)
            @test amplitudes(P) == amp
        end # amplitudes

        @testset "evaluate_lorentzian" begin
            loc = 1.0:10
            amp = 0.1:0.1:0.9
            scl = 1.1
            # single point
            P = PolesContinuedFraction(loc, amp, scl)
            @test evaluate_lorentzian(P, 10, 1) ≈
                0.13282211074263575 - 0.014762307781571633im atol = 10 * eps()
            # grid
            @test evaluate_lorentzian(P, [0.1, 0.3], 0.5) ==
                [evaluate_lorentzian(P, 0.1, 0.5), evaluate_lorentzian(P, 0.3, 0.5)]
        end # evaluate_lorentzian

        @testset "locations" begin
            loc = 0:5
            amp = 6:10
            P = PolesContinuedFraction(loc, amp)
            @test locations(P) == loc
        end # locations

        @testset "scale" begin
            loc = 0:5
            amp = 6:10
            scl = 5
            P = PolesContinuedFraction(loc, amp, scl)
            @test scale(P) === scl
        end # scale

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

    @testset "Core" begin
        @testset "Array" begin
            loc = 1:3
            amp = 4:5
            scl = 2
            P = PolesContinuedFraction(loc, amp, scl)
            @test Array(P) == [1 4 0; 4 2 5; 0 5 3]
            P = PolesContinuedFraction(loc, amp)
            @test Array(P) == [1 4 0; 4 2 5; 0 5 3]
        end # Array
    end # Core

    @testset "Base" begin
        @testset "eltype" begin
            loc = Float64[0.1, 0.2]
            amp = Int[1]
            P = PolesContinuedFraction(loc, amp)
            @test eltype(P) === Float64
        end # eltype

        @testset "length" begin
            loc = [0.1, 0.2]
            amp = [1]
            P = PolesContinuedFraction(loc, amp)
            @test length(P) == 2
        end # length
    end # Base

    @testset "LinearAlgebra" begin
        @testset "SymTridiagonal" begin
            loc = 1:3
            amp = 4:5
            scl = 2
            P = PolesContinuedFraction(loc, amp, scl)
            @test SymTridiagonal(P) == SymTridiagonal(loc, amp)
            P = PolesContinuedFraction(loc, amp)
            @test SymTridiagonal(P) == SymTridiagonal(loc, amp)
        end # SymTridiagonal
    end # LinearAlgebra
end # PolesContinuedFraction
