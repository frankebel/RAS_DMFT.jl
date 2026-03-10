using RAS_DMFT
using LinearAlgebra
using Test

@testset "PolesSumBlock" begin
    @testset "constructor" begin
        loc = 0:1
        wgt = [[1 2; 2 1], [3 4; 4 3]]

        # inner constructor
        P = PolesSumBlock{Int, Int}(loc, wgt)
        @test P.locations == loc
        @test P.weights == wgt
        @test_throws DimensionMismatch PolesSumBlock(rand(3), wgt) # length mismatch
        @test_throws ArgumentError PolesSumBlock([1], [[1 2im; 2im 1]]) # not hermitian
        @test_throws DimensionMismatch PolesSumBlock(0:1, [[1 2im; -2im 1], [1;;]]) # wrong size

        # outer constructors
        P = PolesSumBlock(loc, wgt)
        @test P.locations == loc
        @test P.weights == wgt
        P = PolesSumBlock(loc, [1 + 2im 3im; 4 5 + 6im])
        @test P.locations == loc
        @test P.weights == [[5 4 + 8im; 4 - 8im 16], [9 18 + 15im; 18 - 15im 61]]

        # conversion of type
        P = PolesSumBlock(loc, wgt)
        P_new = PolesSumBlock{UInt, Float64}(P)
        @test typeof(P_new) === PolesSumBlock{UInt, Float64}
        @test P_new.locations == loc
        @test P_new.weights == wgt
    end # constructor

    @testset "custom functions" begin
        @testset "amplitude" begin
            # real
            P = PolesSumBlock(0:1, [[2 0; 0 1], [0 0; 0 0]])
            @inferred amplitude(P, 1)
            @test amplitude(P, 1) == [sqrt(2) 0; 0 1]
            # complex
            P = PolesSumBlock(0:1, [[1 0.5im; -0.5im 1], [0 0; 0 0]])
            @inferred amplitude(P, 1)
            @test norm(
                amplitude(P, 1) -
                    [1 + sqrt(3) (sqrt(3) - 1)im; -(sqrt(3) - 1)im 1 + sqrt(3)] ./ (2 * sqrt(2)),
            ) < 10 * eps()
        end # amplitude

        @testset "amplitudes" begin
            v1 = [1 + 2im, 4] # vector from which first weights are constructed
            v2 = [3im, 5 + 6im] # vector from which second weights are constructed
            P = PolesSumBlock(0:1, [1 + 2im 3im; 4 5 + 6im])
            @inferred amplitudes(P)
            amp = amplitudes(P)
            @test norm(amp[1] - 1 / sqrt(21) * v1 * v1') < 20 * eps()
            @test norm(amp[2] - 1 / sqrt(70) * v2 * v2') < 20 * eps()
        end # amplitudes

        @testset "evaluate_gaussian" begin
            ω = 0.15
            σ = 0.04
            loc = [0.1, 0.2]
            amp = reshape(0.1:0.1:0.4, (2, 2))
            P = PolesSumBlock(loc, amp)
            @test norm(
                evaluate_gaussian(P, ω, σ) - [
                    -1.5277637226549838 - 1.4345225621076145im -1.90970465331873 - 2.00833158695066im
                    -1.90970465331873 - 2.00833158695066im -2.2916455839824765 - 2.8690451242152295im
                ],
            ) < 10 * eps()
        end # evaluate_gaussian

        @testset "evaluate_lorentzian" begin
            ω = 10.0
            δ = 1.0
            loc = 1.0:10
            amp = reshape(0.1:0.1:2, (2, 10))
            # single point
            P = PolesSumBlock(loc, amp)
            @test norm(
                evaluate_lorentzian(P, ω, δ) - [
                    3.419109056634164 - 5.796080126589452im 3.669440321902302 - 6.127487634859227im
                    3.669440321902302 - 6.127487634859227im 3.9413975570979876 - 6.478614061451364im
                ],
            ) < eps()
            # grid
            ω = rand(2)
            @test evaluate_lorentzian(P, ω, δ) ==
                [evaluate_lorentzian(P, ω[1], δ), evaluate_lorentzian(P, ω[2], δ)]
        end # evaluate_lorentzian

        @testset "flip_spectrum!" begin
            P = PolesSumBlock(
                [0.0, 0.25, 0.3],
                [[0.0 0.25; 0.25 0.5], [0.75 1.0; 1.0 2.0], [1.25 1.5; 1.5 3.0]],
            )
            @test flip_spectrum!(P) === P
            @test locations(P) == [-0.3, -0.25, -0.0]
            @test weights(P) ==
                [[1.25 1.5; 1.5 3.0], [0.75 1.0; 1.0 2.0], [0.0 0.25; 0.25 0.5]]
        end # flip_spectrum!

        @testset "merge_degenerate_poles!" begin
            loc = [0.2, 0.3, 0.6]
            wgt = [[1 0; 0 1], [1 0; 0 0], [2 1; 1 2]]
            P = PolesSumBlock(loc, wgt)
            # default tolerance too small
            foo = copy(P)
            @test merge_degenerate_poles!(foo) === foo
            @test locations(foo) == loc
            @test weights(foo) == wgt
            # merge 2 poles
            @test merge_degenerate_poles!(foo, 0.11) === foo
            @test locations(foo) == [0.2, 0.6]
            @test weights(foo) == [[2 0; 0 1], [2 1; 1 2]]
        end # merge_degenerate_poles!

        @testset "merge_negative_locations_to_zero!" begin
            P = PolesSumBlock(
                [-0.1, -0.0, 0.2],
                [[0.0 0.25; 0.25 0.5], [0.75 1.0; 1.0 2.0], [1.25 1.5; 1.5 3.0]],
            )
            @test merge_negative_locations_to_zero!(P) === P
            @test locations(P) == [0.0, 0.2]
            @test weights(P) == [[0.75 1.25; 1.25 2.5], [1.25 1.5; 1.5 3.0]]
        end # merge_negative_locations_to_zero!

        @testset "merge_small_weight!" begin
            loc = [0.0, 1, 4]
            W1 = [2.0 1; 1 3]
            W2 = [4.0 5; 5 6]
            W3 = [7.0 8; 8 9]
            tol = 5.0 # W1 is below, all others above
            # tolerance small
            P = PolesSumBlock([-1.0, 0.0, 1.5], [copy(W1), copy(W2), copy(W3)])
            @test merge_small_weight!(P, eps()) === P
            @test locations(P) == [-1.0, 0.0, 1.5]
            @test weights(P) == [W1, W2, W3]
            # first index
            P = PolesSumBlock(copy(loc), [copy(W1), copy(W2), copy(W3)])
            merge_small_weight!(P, tol)
            @test locations(P) == [1.0, 4.0]
            @test weights(P) == [[6 6; 6 9], [7 8; 8 9]]
            # last index
            P = PolesSumBlock(copy(loc), [copy(W2), copy(W3), copy(W1)])
            merge_small_weight!(P, tol)
            @test locations(P) == [0, 1]
            @test weights(P) == [[4 5; 5 6], [9 9; 9 12]]
            # middle index
            P = PolesSumBlock(copy(loc), [copy(W2), copy(W1), copy(W3)])
            merge_small_weight!(P, tol)
            @test locations(P) == [0, 4]
            @test weights(P) == [[5.5 5.75; 5.75 8.25], [7.5 8.25; 8.25 9.75]]
            # merge zero weight
            P = PolesSumBlock([0, 2], [zeros(2, 2), ones(2, 2)])
            merge_small_weight!(P, 0)
            @test locations(P) == [2]
            @test weights(P) == [ones(2, 2)]
            # Errors
            @test_throws ArgumentError merge_small_weight!(
                PolesSumBlock([4.0, 1.0], [W1, W2]), eps()
            )
        end # merge_small_weight!

        @testset "moment" begin
            P = PolesSumBlock([-0.5, 0.0, 0.5], [0.25 1.5 0.25; 0.5 0.75 2.5])
            @test RAS_DMFT.moment(P) == [2.375 1.875; 1.875 7.0625]
            @test RAS_DMFT.moment(P, 1) == [0.0 0.25; 0.25 3.0]
        end # moment

        @testset "moments" begin
            P = PolesSumBlock([-0.5, 0.0, 0.5], [0.25 1.5 0.25; 0.5 0.75 2.5])
            @test moments(P, 0:1) == [[2.375 1.875; 1.875 7.0625], [0.0 0.25; 0.25 3.0]]
        end # moments

        @testset "remove_zero_weight!" begin
            P = PolesSumBlock([0, 1], [[0 0; 0 0], [1 0; 0 0]])
            # remove zero location
            foo = copy(P)
            @test remove_zero_weight!(foo) === foo
            @test locations(foo) == [1]
            @test weights(foo) == [[1 0; 0 0]]
            # keep zero location
            foo = copy(P)
            @test remove_zero_weight!(foo, false) === foo
            @test locations(foo) == [0, 1]
            @test weights(foo) == [[0 0; 0 0], [1 0; 0 0]]
        end # remove_zero_weight!

        @testset "remove_zero_weight" begin
            P = PolesSumBlock([0, 1], [[0 0; 0 0], [1 0; 0 0]])
            # remove zero location
            foo = remove_zero_weight(P)
            @test locations(foo) == [1]
            @test weights(foo) == [[1 0; 0 0]]
            # keep zero location
            foo = remove_zero_weight(P, false)
            @test locations(foo) == [0, 1]
            @test weights(foo) == [[0 0; 0 0], [1 0; 0 0]]
        end # remove_zero_weight

        @testset "to_grid" begin
            # all poles within grid, middle pole not centered
            P = PolesSumBlock(
                [0.1, 0.25, 0.3],
                [[0.0 0.25; 0.25 0.5], [0.75 1.0; 1.0 2.0], [1.25 1.5; 1.5 3.0]],
            )
            grid = [0.1, 0.3]
            foo = to_grid(P, grid)
            @test locations(foo) == [0.1, 0.3]
            @test norm(weight(foo, 1) - [0.1875 0.5; 0.5 1.0]) < 10 * eps()
            @test norm(weight(foo, 2) - [1.8125 2.25; 2.25 4.5]) < 10 * eps()
        end  # to_grid

        @testset "weight" begin
            loc = 0:1
            wgt = [[1 2; 2 1], [3 4; 4 3]]
            P = PolesSumBlock(loc, wgt)
            @test weight(P, 1) === wgt[1]
            @test weight(P, 2) === wgt[2]
        end # weight

        @testset "weights" begin
            loc = 0:1
            wgt = [[1 2; 2 1], [3 4; 4 3]]
            P = PolesSumBlock(loc, wgt)
            @test weights(P) === wgt
        end # weights
    end # custom functions

    @testset "Core" begin
        @testset "Array" begin
            A = [-7, 9]
            B1 = [4 3; 3 4]
            W1 = B1 * B1'
            B2 = [5 0; 0 1]
            W2 = B2 * B2'
            P = PolesSumBlock(A, [W1, W2])
            m = Array(P)
            @test typeof(m) == Matrix{Float64}
            @test norm(
                m -
                    [
                    0 0  4  3 5 0
                    0 0  3  4 0 1
                    4 3 -7  0 0 0
                    3 4  0 -7 0 0
                    5 0  0  0 9 0
                    0 1  0  0 0 9
                ]
            ) < 50 * eps()
        end # Array
    end # Core

    @testset "Base" begin
        @testset "+" begin
            # addition must sort and merge poles
            A = PolesSumBlock([1, 3], [[1 0; 0 1], [2 1; 1 0]])
            B = PolesSumBlock([3, 2], [[4 5; 5 6], [0 0; 0 0]])
            P = A + B
            @test locations(P) == [1, 2, 3] # must be sorted
            @test weights(P) == [[1 0; 0 1], [0 0; 0 0], [6 6; 6 6]]
            @test weight(P, 1) !== weight(A, 1)
        end

        @testset "copy" begin
            P = PolesSumBlock(rand(2), [rand(1, 1) for _ in 1:2])
            foo = copy(P)
            @test locations(foo) == locations(P)
            @test locations(foo) !== locations(P)
            @test weights(foo) == weights(P)
            @test weights(foo) !== weights(P)
        end # copy

        @testset "eltype" begin
            @test eltype(PolesSumBlock(1:2, [[1 0; 0 0], [0 0; 0 0]])) === Int
            @test eltype(PolesSumBlock(1.0:2, [[1 0; 0 0], [0 0; 0 0]])) === Float64
        end # eltype

        @testset "length" begin
            P = PolesSumBlock(rand(10), rand(4, 10))
            @test length(P) === 10
        end

        @testset "reverse!" begin
            P = PolesSumBlock([1, 2], [[3 4; 4 5], [6 7; 7 8]])
            @test reverse!(P) === P
            @test locations(P) == [2, 1]
            @test weights(P) == [[6 7; 7 8], [3 4; 4 5]]
        end # reverse!

        @testset "reverse" begin
            P = PolesSumBlock([1, 2], [[3 4; 4 5], [6 7; 7 8]])
            foo = reverse(P)
            @test foo !== P
            @test locations(P) == [1, 2]
            @test locations(foo) == [2, 1]
            @test weights(P) == [[3 4; 4 5], [6 7; 7 8]]
            @test weights(foo) == [[6 7; 7 8], [3 4; 4 5]]
        end # reverse

        @testset "size" begin
            P = PolesSumBlock(rand(10), rand(4, 10))
            @test size(P) == (4, 4)
            @test size(P, 1) === 4
            @test size(P, 2) === 4
        end # size

        @testset "sort!" begin
            loc = [2, 1]
            wgt = [[1 0; 0 1], [2 1; 1 0]]
            P = PolesSumBlock(loc, wgt)
            @test sort!(P) === P
            @test locations(P) == [1, 2]
            @test weights(P) == [[2 1; 1 0], [1 0; 0 1]]
        end # sort!

        @testset "sort" begin
            loc = [2, 1]
            wgt = [[1 0; 0 1], [2 1; 1 0]]
            P = PolesSumBlock(loc, wgt)
            foo = sort(P)
            @test foo !== P
            @test locations(P) == [2, 1]
            @test locations(foo) == [1, 2]
            @test weights(P) == [[1 0; 0 1], [2 1; 1 0]]
            @test weights(foo) == [[2 1; 1 0], [1 0; 0 1]]
        end # sort

        @testset "transpose" begin
            P = PolesSumBlock([0, 1], [[5 4 + 8im; 4 - 8im 16], [9 18 + 15im; 18 - 15im 61]])
            Pt = transpose(P)
            @test locations(Pt) == [0, 1]
            @test locations(Pt) !== locations(P) # must copy
            @test weights(Pt) == [[5 4 - 8im; 4 + 8im 16], [9 18 - 15im; 18 + 15im 61]]
        end # transpose
    end # base
end # PolesSumBlock
