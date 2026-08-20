using RAS_DMFT
using LinearAlgebra
using Test

@testset "PolesSumBlock" begin
    @testset "constructors" begin
        @testset "inner constructor" begin
            locs = [0, 1]
            wgts = [[1 2; 2 1], [3 4; 4 3]]
            P = PolesSumBlock{Int, Int}(locs, wgts)
            @test P.locations === locs
            @test P.weights === wgts
            @test @allocated(PolesSumBlock{Int, Int}(locs, wgts)) == 0 # no allocations
            @test_throws DimensionMismatch PolesSumBlock{Int, Int}(rand(Int, 3), wgts) # length mismatch
            @test_throws ArgumentError PolesSumBlock{Int, Complex{Int}}([1], [[1 2im; 2im 1]]) # not Hermitian
            @test_throws DimensionMismatch PolesSumBlock{Int, Int}([0, 1], [[1 2im; -2im 1], [1;;]]) # weights wrong size
            PolesSumBlock{Int, Float64}(Int[], Matrix{Float64}[]) # empty lists
        end # inner constructor

        @testset "outer constructors" begin
            # canonical form
            locs = [0, 1]
            wgts = [[1 2; 2 1], [3 4; 4 3]]
            P = PolesSumBlock(locs, wgts)
            @test P.locations == locs
            @test P.locations !== locs
            @test P.weights == wgts
            @test P.weights !== wgts

            # sort and merge degenerate poles
            locs = [0, 1, 0]
            wgts = [[1 2; 2 1], [3 4; 4 3], [5 6; 6 5]]
            P = PolesSumBlock(locs, wgts)
            @test P.locations == [0, 1]
            @test P.weights == [[6 8; 8 6], [3 4; 4 3]]

            # supply amplitudes
            locs = [0, 1]
            amps = [1 + 2im 3im; 4 5 + 6im]
            P = PolesSumBlock(locs, amps)
            @test P.locations == [0, 1]
            @test P.weights == [[5 4 + 8im; 4 - 8im 16], [9 18 + 15im; 18 - 15im 61]]

            # correct merging of degenerate poles
            locs = [-3.0, -2.5, -2.75]
            amps = [1 5 3; 2 6 4]
            P = PolesSumBlock(locs, amps, 0.25)
            @test P.locations == [-3.0, -2.5]
            @test P.weights == [[10 14; 14 20], [25 30; 30 36]]

            # single pole
            locs = [0.0]
            amps = [1; 2;;]
            P = PolesSumBlock(locs, amps)
            @test P.locations == [0.0]
            @test P.weights == [[1 2; 2 4]]

            # conversion of type
            locs = [0, 1]
            wgts = [[1 2; 2 1], [3 4; 4 3]]
            P = PolesSumBlock(locs, wgts)
            P_new = PolesSumBlock{UInt, Float64}(P)
            @test P_new isa PolesSumBlock{UInt, Float64}
            @test P_new.locations == locs
            @test P_new.weights == wgts
        end # outer constructors
    end # constructors

    @testset "custom functions" begin
        @testset "amplitude" begin
            # real
            P = PolesSumBlock(0:1, [[2.0 0; 0 1], [0 0; 0 0]])
            @inferred amplitude(P, 1)
            @test amplitude(P, 1) == [sqrt(2) 0; 0 1]

            # complex
            P = PolesSumBlock(0:1, [[1 0.5im; -0.5im 1], [0 0; 0 0]])
            @inferred amplitude(P, 1)
            @test norm(
                amplitude(P, 1) -
                    [1 + sqrt(3) (sqrt(3) - 1)im; -(sqrt(3) - 1)im 1 + sqrt(3)] ./ (2 * sqrt(2)),
            ) < 10 * eps()

            # thin with some weights exactly zero
            w = zeros(10, 10)
            w[3, 3] = 4
            w[5, 5] = 9
            amps = zeros(10, 2)
            amps[3, 1] = 2
            amps[5, 2] = 3
            P = PolesSumBlock(rand(1), [w])
            amp2 = amplitude(P, 1, thin = true)
            @test amps == amp2
        end # amplitude

        @testset "amplitudes" begin
            v1 = [1 + 2im, 4] # vector from which first weights are constructed
            v2 = [3im, 5 + 6im] # vector from which second weights are constructed
            P = PolesSumBlock(0:1, [1.0 + 2im 3im; 4 5 + 6im])
            @inferred amplitudes(P)
            amps = amplitudes(P)
            @test norm(amps[1] - 1 / sqrt(21) * v1 * v1') < 20 * eps()
            @test norm(amps[2] - 1 / sqrt(70) * v2 * v2') < 20 * eps()

            # thin with random values
            amps = [rand(ComplexF64, 10, i) for i in 1:10]
            wghts = [amps * amps' for amps in amps]
            locs = sort!(rand(10))
            P = PolesSumBlock(locs, wghts)
            amps2 = amplitudes(P, sqrt(100 * eps()); thin = true)
            for i in eachindex(amps)
                w = wghts[i]
                amps = amps2[i]
                @test size(amps) == (10, i)
                @test norm(w - amps * amps') < 1.0e6 * eps()
            end
        end # amplitudes

        @testset "evaluate" begin
            locs = [-1.0, 0.0, 2.0]
            amps = reshape(0.2:0.1:0.7, (2, 3))
            P = PolesSumBlock(locs, amps)
            # upper/lower complex plane
            z = 0.5 + 1.0im
            @test norm(
                evaluate(P, z) - [
                    -0.083692307692307677 - 0.25107692307692314im -0.08615384615384615 - 0.30769230769230771im
                    -0.08615384615384615 - 0.30769230769230771im -0.084615384615384592 - 0.37846153846153846im
                ],
            ) < 100 * eps()
            @test evaluate(P, conj(z)) == conj.(permutedims(evaluate(P, z)))
            # Matsubara frequency
            z = 2im
            @test norm(
                evaluate(P, z) - [
                    -0.082 - 0.186im -0.093 - 0.229im
                    -0.093 - 0.229im -0.1045 - 0.2835im
                ],
            ) < 100 * eps()
            @test evaluate(P, conj(z)) == conj.(permutedims(evaluate(P, z)))
            # grid
            zs = [0.1 + 0.5im, 0.3 + 0.5im]
            @test evaluate(P, zs) == [evaluate(P, zs[1]), evaluate(P, zs[2])]
        end # evaluate

        @testset "evaluate_gaussian" begin
            ω = 0.5
            σ = 1.0
            locs = [-1.0, 0.0, 2.0]
            amps = reshape(0.2:0.1:0.7, (2, 3))
            P = PolesSumBlock(locs, amps)
            @test norm(
                evaluate_gaussian(P, ω, σ) - [
                    -0.16702850198727615 - 0.33972394588525767im -0.178700179083296 - 0.41651710181548046im
                    -0.178700179083296 - 0.41651710181548046im -0.18576841335312091 - 0.51250854672825896im
                ],
            ) < 10 * eps()
            # grid
            ω = [0.1, 0.3]
            @test evaluate_gaussian(P, ω, 0.5) ==
                [evaluate_gaussian(P, ω[1], 0.5), evaluate_gaussian(P, ω[2], 0.5)]
        end # evaluate_gaussian

        @testset "evaluate_lorentzian" begin
            ω = 0.5
            δ = 1.0
            locs = [-1.0, 0.0, 2.0]
            amps = reshape(0.2:0.1:0.7, (2, 3))
            # single point
            P = PolesSumBlock(locs, amps)
            @test norm(
                evaluate_lorentzian(P, ω, δ) - [
                    -0.083692307692307677 - 0.25107692307692314im -0.08615384615384615 - 0.30769230769230771im
                    -0.08615384615384615 - 0.30769230769230771im -0.084615384615384592 - 0.37846153846153846im
                ],
            ) < eps()
            # grid
            ω = [0.1, 0.3]
            @test evaluate_lorentzian(P, ω, 0.5) ==
                [evaluate_lorentzian(P, ω[1], 0.5), evaluate_lorentzian(P, ω[2], 0.5)]
        end # evaluate_lorentzian

        @testset "filling" begin
            P = PolesSumBlock(
                -1:1,
                [[1 2; 2 3], [4 5; 5 6], [7 8; 8 9]],
            )
            @inferred filling(P)
            @test filling(P) == [3 4.5; 4.5 6]
            @test filling(P, -Inf) == zeros(2, 2)
            @test filling(P, -1.1) == zeros(2, 2)
            @test filling(P, -1) == [0.5 1; 1 1.5]
            @test filling(P, 0) == [3 4.5; 4.5 6]
            @test filling(P, 1) == [8.5 11; 11 13.5]
            @test filling(P, 1.1) == [12 15; 15 18]
            @test filling(P, Inf) == [12 15; 15 18]
        end # filling

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

        @testset "location" begin
            P = PolesSumBlock(
                [0.0, 0.25, 0.3],
                [[0.0 0.25; 0.25 0.5], [0.75 1.0; 1.0 2.0], [1.25 1.5; 1.5 3.0]],
            )
            @test location(P, 1) == 0.0
            @test location(P, 2) == 0.25
            @test location(P, 3) == 0.3
            @test_throws BoundsError location(P, 4)
        end # location

        @testset "locations" begin
            P = PolesSumBlock(
                [0.0, 0.25, 0.3],
                [[0.0 0.25; 0.25 0.5], [0.75 1.0; 1.0 2.0], [1.25 1.5; 1.5 3.0]],
            )
            @test locations(P) == [0.0, 0.25, 0.3]
        end # locations

        @testset "merge_degenerate_poles!" begin
            locs = [0.2, 0.3, 0.6]
            wgts = [[1 0; 0 1], [1 0; 0 0], [2 1; 1 2]]
            P = PolesSumBlock(locs, wgts)
            # default tolerance too small
            foo = copy(P)
            @test merge_degenerate_poles!(foo) === foo
            @test locations(foo) == locs
            @test weights(foo) == wgts
            # merge 2 poles
            @test merge_degenerate_poles!(foo, 0.11) === foo
            @test locations(foo) == [0.2, 0.6]
            @test weights(foo) == [[2 0; 0 1], [2 1; 1 2]]
            # merge poles within [-tol, tol] into zero
            P = PolesSumBlock([-0.05, 0.05, 0.6], [[1 0; 0 1], [1 0; 0 0], [2 1; 1 2]])
            @test merge_degenerate_poles!(P, 0.1) === P
            @test locations(P) == [0.0, 0.6]
            @test weights(P) == [[2 0; 0 1], [2 1; 1 2]]
            # merge negative poles
            P = PolesSumBlock([-0.8, -0.6, 0.4], [[1 0; 0 1], [1 0; 0 0], [2 1; 1 2]])
            @test merge_degenerate_poles!(P, 0.25) === P
            @test locations(P) == [-0.6, 0.4]
            @test weights(P) == [[2 0; 0 1], [2 1; 1 2]]
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
            locs = [0.0, 1, 4]
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
            P = PolesSumBlock(copy(locs), [copy(W1), copy(W2), copy(W3)])
            merge_small_weight!(P, tol)
            @test locations(P) == [1.0, 4.0]
            @test weights(P) == [[6 6; 6 9], [7 8; 8 9]]
            # last index
            P = PolesSumBlock(copy(locs), [copy(W2), copy(W3), copy(W1)])
            merge_small_weight!(P, tol)
            @test locations(P) == [0, 1]
            @test weights(P) == [[4 5; 5 6], [9 9; 9 12]]
            # middle index
            P = PolesSumBlock(copy(locs), [copy(W2), copy(W1), copy(W3)])
            merge_small_weight!(P, tol)
            @test locations(P) == [0, 4]
            @test weights(P) == [[5.5 5.75; 5.75 8.25], [7.5 8.25; 8.25 9.75]]
            # merge zero weight
            P = PolesSumBlock([0, 2], [zeros(2, 2), ones(2, 2)])
            merge_small_weight!(P, 0)
            @test locations(P) == [2]
            @test weights(P) == [ones(2, 2)]
        end # merge_small_weight!

        @testset "moment" begin
            P = PolesSumBlock([-0.5, 0.0, 0.5], [0.25 1.5 0.25; 0.5 0.75 2.5])
            @test moment(P) == [2.375 1.875; 1.875 7.0625]
            @test moment(P, 1) == [0.0 0.25; 0.25 3.0]
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

        @testset "shift_spectrum!" begin
            P = PolesSumBlock([0, 1], [[0 0; 0 0], [1 0; 0 0]])
            @test RAS_DMFT.shift_spectrum!(P, 2) === P
            @test locations(P) == [-2, -1]
            @test weights(P) == [[0 0; 0 0], [1 0; 0 0]]
        end # shift_spectrum!

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
            locs = 0:1
            wgts = [[1 2; 2 1], [3 4; 4 3]]
            P = PolesSumBlock(locs, wgts)
            @test weight(P, 1) == wgts[1]
            @test weight(P, 2) == wgts[2]
        end # weight

        @testset "weights" begin
            locs = 0:1
            wgts = [[1 2; 2 1], [3 4; 4 3]]
            P = PolesSumBlock(locs, wgts)
            @test weights(P) == wgts
        end # weights
    end # custom functions

    @testset "Base" begin
        @testset "+" begin
            # addition must merge degenerate poles
            A = PolesSumBlock([1, 3], [[1 0; 0 1], [2 1; 1 0]])
            B = PolesSumBlock([2, 3], [[0 0; 0 0], [4 5; 5 6]])
            P = A + B
            @test locations(P) == [1, 2, 3]
            @test weights(P) == [[1 0; 0 1], [0 0; 0 0], [6 6; 6 6]]
            @test weight(P, 1) !== weight(A, 1)
            @test weight(P, 2) !== weight(B, 1)

            # B exhausted before A -> append remaining A-poles
            A = PolesSumBlock([1, 3], [[1 0; 0 1], [2 1; 1 0]])
            B = PolesSumBlock([2], [[4 5; 5 6]])
            P = A + B
            @test locations(P) == [1, 2, 3]
            @test weights(P) == [[1 0; 0 1], [4 5; 5 6], [2 1; 1 0]]
        end

        @testset "convert" begin
            P = PolesSumBlock([1, 3], [[1 0; 0 1], [2 1; 1 0]])
            P_new = convert(PolesSumBlock{Float64, ComplexF64}, P)
            @test P_new isa PolesSumBlock{Float64, ComplexF64}
            @test locations(P_new) == locations(P)
            @test weights(P_new) == weights(P)
            P_new = convert(PolesSumBlock{Int, Int}, P)
            @test P_new === P
        end # convert

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

        @testset "isempty" begin
            @test isempty(PolesSumBlock(Int[], Matrix{Float64}[]))
            @test !isempty(PolesSumBlock(rand(10), rand(2, 10)))
        end # isempty

        @testset "iterate" begin
            P = PolesSumBlock(1:2, [[1 0; 0 0], [0 0; 0 0]])
            @test iterate(P) == ((1, [1 0; 0 0]), 1)
            @test iterate(P, 1) == ((2, [0 0; 0 0]), 2)
            @test iterate(P, 2) === nothing
        end # iterate

        @testset "length" begin
            @test length(PolesSumBlock(Int[], Matrix{Float64}[])) === 0
            @test length(PolesSumBlock(rand(10), rand(4, 10))) === 10
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

        @testset "show" begin
            P = PolesSumBlock(Int[], Matrix{Float64}[])
            @test sprint(show, P) == "PolesSumBlock{Int64, Float64} with 0 poles"
            P = PolesSumBlock(rand(Int, 1), [hermitianpart!(rand(Float64, 2, 2))])
            @test sprint(show, P) == "PolesSumBlock{Int64, Float64} with 1 poles of size 2×2"
            P = PolesSumBlock(rand(Int, 2), [hermitianpart!(rand(Float64, 3, 3)) for _ in 1:2])
            @test sprint(show, P) == "PolesSumBlock{Int64, Float64} with 2 poles of size 3×3"
        end # show

        @testset "size" begin
            P = PolesSumBlock(rand(10), rand(4, 10))
            @test size(P) == (4, 4)
            @test size(P, 1) === 4
            @test size(P, 2) === 4
        end # size

        @testset "sort!" begin
            locs = [2, 1]
            wgts = [[1 0; 0 1], [2 1; 1 0]]
            P = PolesSumBlock(locs, wgts)
            @test locations(P) == [1, 2]
            @test weights(P) == [[2 1; 1 0], [1 0; 0 1]]
        end # sort!

        @testset "transpose" begin
            P = PolesSumBlock([0, 1], [[5 4 + 8im; 4 - 8im 16], [9 18 + 15im; 18 - 15im 61]])
            Pt = transpose(P)
            @test locations(Pt) == [0, 1]
            @test locations(Pt) !== locations(P) # must copy
            @test weights(Pt) == [[5 4 - 8im; 4 + 8im 16], [9 18 - 15im; 18 + 15im 61]]
        end # transpose
    end # base

    @testset "LinearAlgebra" begin
        @testset "rmul!" begin
            P = PolesSumBlock([-1.0, 0.0], [[1 2; 2 4], [4 -3; -3 7]])
            @test rmul!(P, 2) === P
            @test locations(P) == [-1.0, 0.0] # unchanged
            @test weights(P) == [[2 4; 4 8], [8 -6; -6 14]]
        end # rmul!

        @testset "tr" begin
            P = PolesSumBlock([-1, 0], [[1 + 0.0im 2; 2 4], [4 -3; -3 7]])
            Ps = tr(P)
            @test Ps isa PolesSum{Int, Float64}
            @test locations(Ps) == locations(P)
            @test locations(Ps) !== locations(P)
            @test weights(Ps) == [5.0, 11.0]
        end # tr
    end # LinearAlgebra
end # PolesSumBlock
