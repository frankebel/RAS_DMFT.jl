using RAS_DMFT
using Distributions
using LinearAlgebra
using Test

@testset "PolesSum" begin
    @testset "constructor" begin
        loc = 0:5
        wgt = 5:10

        # inner constructor
        P = PolesSum{Int, Int}(loc, wgt)
        @test typeof(P) === PolesSum{Int, Int}
        @test P.locations == loc
        @test P.weights == wgt
        # length mismatch
        @test_throws DimensionMismatch PolesSum{Int, Int}(rand(3), rand(4))
        @test_throws DimensionMismatch PolesSum{Int, Int}(rand(4), rand(3))

        # outer constructor
        P = PolesSum(loc, wgt)
        @test P.locations == loc
        @test P.weights == wgt

        # conversion of type
        P = PolesSum(loc, wgt)
        P_new = PolesSum{UInt, Float64}(P)
        @test typeof(P_new) === PolesSum{UInt, Float64}
        @test P_new.locations == loc
        @test P_new.weights == wgt
    end # constructor

    @testset "custom functions" begin
        @testset "add_pole_at_zero!" begin
            # pole already exists
            P = PolesSum([-1, 0, 2], [3, 4, 5])
            @test add_pole_at_zero!(P) === P
            @test locations(P) == [-1, 0, 2]
            @test weights(P) == [3, 4, 5]
            # pole does not exist
            P = PolesSum([-1, 2], [3, 5])
            @test add_pole_at_zero!(P) === P
            @test locations(P) == [-1, 0, 2]
            @test weights(P) == [3, 0, 5]
        end # add_pole_at_zero!

        @testset "amplitude" begin
            loc = 0:5
            wgt = 5:10
            P = PolesSum(loc, wgt)
            @test_throws BoundsError amplitude(P, 0)
            @test amplitude(P, 1) == sqrt(5)
            @test amplitude(P, 2) == sqrt(6)
            @test amplitude(P, 3) == sqrt(7)
            @test amplitude(P, 4) == sqrt(8)
            @test amplitude(P, 5) == sqrt(9)
            @test amplitude(P, 6) == sqrt(10)
            @test_throws BoundsError amplitude(P, 7)
        end # amplitude

        @testset "amplitudes" begin
            loc = 0:5
            wgt = 5:10
            P = PolesSum(loc, wgt)
            @test amplitudes(P) == sqrt.(5:10)
            # errors
            P = PolesSum(loc, rand(ComplexF64, 6))
            @test_throws MethodError amplitudes(P)
            P = PolesSum(loc, -wgt)
            @test_throws DomainError amplitudes(P)
        end # amplitudes

        @testset "evaluate_gaussian" begin
            loc = [0.1, 0.2]
            wgt = [0.01, 0.09]
            P = PolesSum(loc, wgt)
            # single point
            ω = 0.15
            σ = 0.04
            @test evaluate_gaussian(P, ω, σ) ≈ -1.5277637226549838 - 1.4345225621076145im atol =
                10 * eps()
            # semicircular DOS
            G = greens_function_bethe_simple(3001)
            ω = -3:0.01:3
            σ = 0.01
            # constant broadening
            h = evaluate_gaussian(G, ω, σ)
            ex = π .* pdf.(Semicircle(1), ω) # exact solution
            @test norm(ex + imag(h)) < 0.2
            @test maximum(abs.(ex + imag(h))) < 0.12
            @test findmin(imag(h))[2] == cld(length(ω), 2) # symmetric
        end # evaluate_gaussian

        @testset "evaluate_lorentzian" begin
            loc = 1.0:10
            wgt = abs2.(0.1:0.1:1)
            P = PolesSum(loc, wgt)
            # single point
            @test evaluate_lorentzian(P, 10, 1) ≈ 0.9853493892744969 - 1.619653515362841im atol =
                10 * eps()
            # grid
            @test evaluate_lorentzian(P, [0.1, 0.3], 0.5) ==
                [evaluate_lorentzian(P, 0.1, 0.5), evaluate_lorentzian(P, 0.3, 0.5)]
        end # evaluate_lorentzian

        @testset "flip_spectrum!" begin
            P = PolesSum([0.1, 0.2], [0.3, 0.4])
            @test flip_spectrum!(P) === P
            @test locations(P) == [-0.2, -0.1]
            @test weights(P) == [0.4, 0.3]
        end # flip_spectrum!

        @testset "flip_spectrum" begin
            P = PolesSum([0.1, 0.2], [0.3, 0.4])
            foo = flip_spectrum(P)
            @test foo !== P
            @test locations(foo) == [-0.2, -0.1]
            @test weights(foo) == [0.4, 0.3]
        end # flip_spectrum

        @testset "locations" begin
            P = PolesSum(0:5, 5:10)
            @test locations(P) === P.locations
        end # locations

        @testset "merge_degenerate_poles!" begin
            P = PolesSum([0.2, 0.3, 0.6], [0.0625, 0.5625, 2.25])
            # manual tolerance
            P1 = copy(P)
            @test merge_degenerate_poles!(P1, 0.11) === P1
            @test locations(P1) == [0.2, 0.6]
            @test weights(P1) == [0.625, 2.25]
            # default tolerance too small
            P1 = copy(P)
            merge_degenerate_poles!(P1)
            @test locations(P1) == [0.2, 0.3, 0.6]
            @test weights(P1) == [0.0625, 0.5625, 2.25]
            # custom tolerance
            P1 = copy(P)
            locations(P1)[2] = 0.5999999999999
            merge_degenerate_poles!(P1, 1.0e-10)
            @test locations(P1) == [0.2, 0.5999999999999]
            @test weights(P1) == [0.0625, 2.8125]
            # negative locations
            P = PolesSum([-0.6, -0.3, -0.2], [2.25, 0.5625, 0.0625])
            merge_degenerate_poles!(P, 0.11)
            @test locations(P) == [-0.6, -0.2]
            @test weights(P) == [2.25, 0.625]
            # poles around zero
            P = PolesSum([-0.03, -0.01, 0.01, 0.05], [0.0625, 0.5625, 2.25, 6.25])
            merge_degenerate_poles!(P, 0.04)
            @test locations(P) == [0.0, 0.05]
            @test weights(P) == [2.875, 6.25]
            # poles at exactly same location
            P = PolesSum([-0.5, -0.5, 0.0, 0.05], [0.0625, 0.5625, 2.25, 6.25])
            merge_degenerate_poles!(P)
            @test locations(P) == [-0.5, 0.0, 0.05]
            @test weights(P) == [0.625, 2.25, 6.25]
            # Errors
            @test_throws ArgumentError merge_degenerate_poles!(
                PolesSum([0.0, -0.1, 0.5], rand(3))
            )
        end # merge_degenerate_poles!

        @testset "merge_negative_locations_to_zero!" begin
            P = PolesSum([-0.1, -0.0, 0.0, 0.2], [0.25, 0.5625, 6.25, 1.5])
            @test merge_negative_locations_to_zero!(P) === P
            @test locations(P) == [0.0, 0.2]
            @test weights(P) == [7.0625, 1.5]
            # degeneracy at zero
            P = PolesSum([-0.0, -0.0, 0.0, 0.2], [0.25, 0.5625, 6.25, 1.5])
            @test merge_negative_locations_to_zero!(P) === P
            @test locations(P) == [0.0, 0.2]
            @test weights(P) == [7.0625, 1.5]
            # # no negative location
            P = PolesSum([0.1, 0.1, 0.5, 1.0], [0.5, 0.75, 2.5, 1.5])
            @test merge_negative_locations_to_zero!(P) === P
            @test locations(P) == [0.1, 0.1, 0.5, 1.0]
            @test weights(P) == [0.5, 0.75, 2.5, 1.5]
            # Errors
            @test_throws ArgumentError merge_negative_locations_to_zero!(
                PolesSum([0.0, -0.1, 0.5], rand(3))
            )
        end # merge_negative_locations_to_zero!

        @testset "merge_negative_weight!" begin
            # equidistant grid
            loc = [-0.5, 0.5, 1.5]
            wgt = [1.5, -0.5, 5.0]
            P = merge_negative_weight!(PolesSum(loc, wgt))
            @test locations(P) === loc
            @test weights(P) === wgt
            @test loc == [-0.5, 0.5, 1.5]
            @test wgt == [1.25, 0.0, 4.75]
            # not equidistant grid
            loc = [-0.5, 0.0, 1.5]
            wgt = [1.5, -0.5, 5.0]
            merge_negative_weight!(PolesSum(loc, wgt))
            @test loc == [-0.5, 0.0, 1.5]
            @test wgt == [1.125, 0.0, 4.875]
            # first pole negative
            loc = [0.0, 1.0, 5.0]
            wgt = [-1.0, 0.5, 2.25]
            merge_negative_weight!(PolesSum(loc, wgt))
            @test loc == [0.0, 1.0, 5.0]
            @test wgt == [0.0, 0.0, 1.75]
            # last pole negative
            loc = [0.0, 1.0, 5.0]
            wgt = [2.25, 0.5, -1.0]
            merge_negative_weight!(PolesSum(loc, wgt))
            @test loc == [0.0, 1.0, 5.0]
            @test wgt == [1.75, 0.0, 0.0]
            # weight exactly cancel
            loc = [0.0, 1.0, 5.0]
            wgt = [-1.0, 0.5, 0.5]
            merge_negative_weight!(PolesSum(loc, wgt))
            @test loc == [0.0, 1.0, 5.0]
            @test wgt == [0.0, 0.0, 0.0]
            # symmetric case
            loc = [-2.0, -0.5, 0.0, 0.5, 2.0]
            wgt = [5.0, -2.0, 1.0, -2.0, 5.0]
            merge_negative_weight!(PolesSum(loc, wgt))
            @test loc == [-2.0, -0.5, 0.0, 0.5, 2.0]
            @test norm(wgt - [3.5, 0.0, 0.0, 0.0, 3.5]) < 10 * eps()
            # previous pole would get negative weight
            loc = [-1.0, -0.5, 0.0, 1.5]
            wgt = [2.0, 1.5, -2.5, 5.0]
            merge_negative_weight!(PolesSum(loc, wgt))
            @test loc == [-1.0, -0.5, 0.0, 1.5]
            @test wgt == [1.7, 0.0, 0.0, 4.3]
        end # merge_negative_weight!

        @testset "merge_small_weight!" begin
            P = PolesSum([-1.0, 0.0, 1.5], [0.25, 6.25, 5.0])
            # tolerance small
            @test merge_small_weight!(P, eps()) === P
            @test locations(P) == [-1.0, 0.0, 1.5]
            @test weights(P) == [0.25, 6.25, 5.0]
            # first index
            merge_small_weight!(P, 1.0)
            @test locations(P) == [0.0, 1.5]
            @test weights(P) == [6.5, 5.0]
            # last index
            P = PolesSum([-1.0, 0.0, 1.5], [5.0, 6.25, 0.25])
            merge_small_weight!(P, 1.0)
            @test locations(P) == [-1.0, 0.0]
            @test weights(P) == [5.0, 6.5]
            # middle index
            P = PolesSum([-1.0, 0.0, 1.5], [25.0, 0.25, 6.25])
            merge_small_weight!(P, 1.0)
            @test locations(P) == [-1.0, 1.5]
            @test weights(P) == [25.15, 6.35]
            # merge zero weight
            P = PolesSum([0, 2], [0, 1])
            merge_small_weight!(P, 0)
            @test locations(P) == [2]
            @test weights(P) == [1]
            # Errors
            @test_throws ArgumentError merge_small_weight!(
                PolesSum([0.0, -0.1, 0.5], rand(3)), eps()
            )
        end # merge_small_weight!

        @testset "moment" begin
            P = PolesSum([-0.5, 0.0, 0.5], [0.0625, 2.25, 0.0625])
            @test RAS_DMFT.moment(P) == 2.375
            @test iszero(RAS_DMFT.moment(P, 1))
            @test RAS_DMFT.moment(P, 2) == 0.03125
            @test iszero(RAS_DMFT.moment(P, 101))
            # odd moment must vanish for even function
            P = PolesSum([-1.0, -eps(), -2.0, 2.0, eps(), 1.0], fill(1.0, 6))
            @test iszero(RAS_DMFT.moment(P, 1))
            @test iszero(RAS_DMFT.moment(P, 101))
        end # moment

        @testset "moments" begin
            P = PolesSum([-0.5, 0.0, 0.5], [0.0625, 2.25, 0.0625])
            @test moments(P, 0:1) == [2.375, 0]
            @test all(iszero, moments(P, 1:2:11))
        end # moments

        @testset "remove_zero_weight!" begin
            P = PolesSum(1:6, [0, 7, 0, 9, 0, -0.0])
            @test remove_zero_weight!(P) === P
            @test locations(P) == [2, 4]
            @test weights(P) == [7, 9]
            # pole at origin
            loc = [-1, 0, 1]
            wgt = [2, 0, 0]
            P = PolesSum(copy(loc), copy(wgt))
            remove_zero_weight!(P)
            @test locations(P) == [-1]
            @test weights(P) == [2]
            P = PolesSum(copy(loc), copy(wgt))
            remove_zero_weight!(P, false)
            @test locations(P) == [-1, 0]
            @test weights(P) == [2, 0]
        end # remove_zero_weight!

        @testset "remove_zero_weight" begin
            P = PolesSum(1:6, [0, 7, 0, 9, 0, -0.0])
            P_new = remove_zero_weight(P)
            @test P_new !== P
            @test locations(P_new) == [2, 4]
            @test weights(P_new) == [7, 9]
            @test locations(P) == 1:6
            @test weights(P) == [0, 7, 0, 9, 0, 0]
        end # remove_zero_weight

        @testset "spectral_function_loggauss" begin
            Λ = 1.2
            N = 150
            loc = grid_log(1, Λ, N)
            grid = [-reverse(loc); 0; loc]
            G = greens_function_bethe_grid(grid)
            W = range(-2; stop = 2, length = 4000) # exclude ω == 0
            P = spectral_function_loggaussian(G, W, 0.2)
            P .*= π
            @test all(>=(0), P) # positive semidefinite
            @test norm(P - reverse(P)) / 4000 < 10 * eps() # symmetry
            @test abs(P[2000] - 2) < 0.01 # Luttinger pinning
            @test first(P) < 1.0e-6 # decay for ω → ±∞
            @test last(P) < 1.0e-6 # decay for ω → ±∞
        end # spectral_function_loggauss

        @testset "to_grid" begin
            # all poles within grid, middle pole centered
            P = PolesSum([0.1, 0.2, 0.3], [25.0, 100.0, 1.0])
            grid = [0.1, 0.3]
            foo = to_grid(P, grid)
            @test locations(foo) !== grid
            @test locations(foo) == [0.1, 0.3]
            @test norm(weights(foo) - [75, 51]) < 100 * eps()
            # all poles within grid, middle pole not centered
            P = PolesSum([0.1, 0.25, 0.3], [25.0, 100.0, 1.0])
            grid = [0.1, 0.3]
            foo = to_grid(P, grid)
            @test locations(foo) == [0.1, 0.3]
            @test weights(foo) == [50, 76]
            # pole outside grid
            P = PolesSum([0.0, 1.0], [25.0, 100.0])
            grid = [0.1, 0.3]
            foo = to_grid(P, grid)
            @test locations(foo) == [0.1, 0.3]
            @test weights(foo) == [25.0, 100.0]
            # poles very close to grid
            P = PolesSum([4.0e-16, 0.9999999999999998], [16.0, 25.0])
            grid = [0.0, 1.0]
            foo = to_grid(P, grid)
            @test locations(foo) == [0.0, 1.0]
            @test weights(foo) == [16.0, 25.0]
        end  # to_grid

        @testset "weight" begin
            P = PolesSum([-1.0, 0.0, 0.5], [0.25, 1.5, 2.5])
            @test_throws BoundsError weight(P, 0)
            @test weight(P, 1) == 0.25
            @test weight(P, 2) == 1.5
            @test weight(P, 3) == 2.5
            @test_throws BoundsError weight(P, 4)
        end # weight

        @testset "weights" begin
            P = PolesSum(0:5, 5:10)
            @test weights(P) === P.weights
        end # weights
    end # custom functions

    @testset "Core" begin
        @testset "Array" begin
            loc = 1:5
            amp = 6:10
            P = PolesSum(loc, abs2.(amp))
            m = Array(P)
            @test typeof(m) === Matrix{Int}
            @test m == [
                0 6 7 8 9 10
                6 1 0 0 0 0
                7 0 2 0 0 0
                8 0 0 3 0 0
                9 0 0 0 4 0
                10 0 0 0 0 5
            ]
            # poles with zero weight
            loc = 1:5
            amp = [6, 7, 0, 9, 0]
            P = PolesSum(loc, abs2.(amp))
            m = Array(P)
            @test m == [
                0 6 7 0 9 0
                6 1 0 0 0 0
                7 0 2 0 0 0
                0 0 0 3 0 0
                9 0 0 0 4 0
                0 0 0 0 0 5
            ]
            # correct promotion
            loc = 1:2
            amp = [1.1, 5.5]
            P = PolesSum(loc, abs2.(amp))
            m = Array(P)
            @test m == [
                0 1.1 5.5
                1.1 1 0
                5.5 0 2
            ]
        end # Array
    end # Core

    @testset "Base" begin
        @testset "+" begin
            # addition must sort resulting poles
            A = PolesSum([0.1, 0.3, 0.2], [0.1, 0.3, 0.2])
            B = PolesSum([0.6, 0.5, 0.4], [6, 5, 4])
            P = A + B
            @test locations(P) == [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
            @test weights(P) == [0.1, 0.2, 0.3, 4.0, 5.0, 6.0]
            # addition must merge degenerate poles
            A = PolesSum([0.1, 0.2], [0.1, 0.25])
            B = PolesSum([0.2], [1])
            P = A + B
            @test locations(P) == [0.1, 0.2]
            @test weights(P) == [0.1, 1.25]
        end

        @testset "-" begin
            # subtraction must sort resulting poles
            A = PolesSum([0.1, 0.3, 0.2], [0.1, 0.3, 0.2])
            B = PolesSum([0.6, 0.5, 0.4], [6, 5, 4])
            P = A - B
            @test locations(P) == [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
            @test weights(P) == [0.1, 0.2, 0.3, -4.0, -5.0, -6.0]
            # addition must merge degenerate poles
            A = PolesSum([0.1, 0.2], [0.1, 0.25])
            B = PolesSum([0.2], [1])
            P = A - B
            @test locations(P) == [0.1, 0.2]
            @test weights(P) == [0.1, -0.75]
        end # -

        @testset "allunique" begin
            @test allunique(PolesSum([0.1, 0.0, -0.5], rand(3)))
            @test !allunique(PolesSum([0.1, 0.0, 0.1], rand(3)))
            @test !allunique(PolesSum([-0.1, 0.0, -0.1], rand(3)))
            @test !allunique(PolesSum([0.0, 0.0, -0.1], rand(3)))
            @test !allunique(PolesSum([-0.0, 0.0, -0.1], rand(3)))
        end # allunique

        @testset "copy" begin
            loc = 1:5
            wgt = 6:10
            P = PolesSum(loc, wgt)
            foo = copy(P)
            @test typeof(foo) === typeof(P)
            @test locations(foo) !== locations(P)
            @test locations(foo) == locations(P)
            @test weights(foo) !== weights(P)
            @test weights(foo) == weights(P)
        end # copy

        @testset "eltype" begin
            @test eltype(PolesSum([0, 1], [0, 1])) === Int64
            @test eltype(PolesSum([0, 1], [0.0, 1.0])) === Float64
            @test eltype(PolesSum([0.0, 1.0], [0.0im, 1.0])) === ComplexF64
        end # eltype

        @testset "inv" begin
            grid = range(-1, 1; length = 101)
            G = greens_function_bethe_grid(grid)
            a0, P = inv(G)
            @test length(P) === 100 # originally 101 poles
            # poles are symmetric
            @test abs(a0) < eps()
            @test norm(locations(P) + reverse(locations(P))) < 50 * eps()
            @test norm(weights(P) - reverse(weights(P))) < 10 * eps()
            @test RAS_DMFT.moment(P, 0) ≈ 0.25 atol = 1.0e-4 # total weight
            # evaluate
            δ = 0.1
            @test norm(
                evaluate_lorentzian(G, 0, δ) -
                    1 / (im * δ - a0 - evaluate_lorentzian(P, 0, δ)),
            ) < 50 * eps()
            ω = 1.0
            δ = 0.1
            @test norm(
                evaluate_lorentzian(G, ω, δ) -
                    1 / (ω + im * δ - a0 - evaluate_lorentzian(P, ω, δ)),
            ) < 10 * eps()
            # symmetry
            z1 = 1 / (-0.8 + 0.1im - a0 - evaluate_lorentzian(P, -0.8, 0.1))
            z2 = 1 / (0.8 + 0.1im - a0 - evaluate_lorentzian(P, 0.8, 0.1))
            @test real(z1) ≈ -real(z2) rtol = 20 * eps()
            @test imag(z1) ≈ imag(z2) rtol = 20 * eps()
        end # inv

        @testset "issorted" begin
            @test issorted(PolesSum([-0.3, 0.0, 0.1], rand(3)))
            @test issorted(PolesSum([-0.0, 0.0, 0.1], rand(3)))
            @test issorted(PolesSum([0.0, 0.0, 0.1], rand(3)))
            @test !issorted(PolesSum([0.0, -0.0, 0.1], rand(3)))
            @test !issorted(PolesSum([0.0, 0.2, 0.1], rand(3)))
            @test issorted(PolesSum([0.2, 0.1, 0.0], rand(3)); rev = true)
        end # issorted

        @testset "length" begin
            @test length(PolesSum(rand(10), rand(10))) === 10
        end # length

        @testset "reverse!" begin
            P = PolesSum([0.1, 0.2], [0.3, 0.4])
            @test reverse!(P) === P
            @test locations(P) == [0.2, 0.1]
            @test weights(P) == [0.4, 0.3]
        end # reverse!

        @testset "reverse" begin
            P = PolesSum([0.1, 0.2], [0.3, 0.4])
            foo = reverse(P)
            @test foo !== P
            @test locations(P) == [0.1, 0.2]
            @test locations(foo) == [0.2, 0.1]
            @test weights(P) == [0.3, 0.4]
            @test weights(foo) == [0.4, 0.3]
        end # reverse

        @testset "sort!" begin
            loc = [2, 1]
            wgt = [9, 16]
            P = PolesSum(loc, wgt)
            @test sort!(P) === P
            @test locations(P) == [1, 2]
            @test weights(P) == [16, 9]
        end # sort!

        @testset "sort" begin
            loc = [2, 1]
            wgt = [3, 4]
            P = PolesSum(loc, wgt)
            foo = sort(P)
            @test foo !== P
            @test locations(P) == [2, 1]
            @test locations(foo) == [1, 2]
            @test weights(P) == [3, 4]
            @test weights(foo) == [4, 3]
        end # sort
    end # Base

    @testset "LinearAlgebra" begin
        @testset "axpby!" begin
            x = PolesSum([-1.0, 0.0, 1.0], [0.5, 0.75, 2.0])
            y = PolesSum([-2.0, 0.0, 2.0], [0.75, 2.0, 1.0])
            @test axpby!(0.5, x, 2.0, y) === y
            @test locations(y) == [-2.0, -1.0, 0.0, 1.0, 2.0]
            @test weights(y) == [1.5, 0.25, 4.375, 1.0, 2.0]
            # x must be unchanged
            @test locations(x) == [-1.0, 0.0, 1.0]
            @test weights(x) == [0.5, 0.75, 2.0]
        end # axpby!

        @testset "rmul!" begin
            P = PolesSum([-1.0, 0.0, 2.0], [0.5, 1.2, 2.0])
            @test rmul!(P, 2) === P
            @test locations(P) == [-1.0, 0.0, 2.0] # unchanged
            @test weights(P) == [1.0, 2.4, 4.0]
        end # rmul!
    end # LinearAlgebra
end # PolesSum
