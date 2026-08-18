using RAS_DMFT
using LinearAlgebra
using Test

@testset "Green's function" begin
    @testset "Bethe lattice" begin
        @testset "analytic" begin
            # test different number types
            @test greens_function_bethe_analytic(-2) == -0.5358983848622456
            @test greens_function_bethe_analytic(-true) == -2
            @test greens_function_bethe_analytic(-1 // 2) == -1.0 - 1.7320508075688772im
            @test greens_function_bethe_analytic(false) == -2im
            @test greens_function_bethe_analytic(0.1im) == -1.809975124224178im
            @test greens_function_bethe_analytic(0.5) == 1.0 - 1.7320508075688772im
            @test greens_function_bethe_analytic(0x01) == 2
            @test greens_function_bethe_analytic(2.0) == 0.5358983848622456
            # vary half-bandwidth
            @test greens_function_bethe_analytic(-1, 2) == -0.5 - 0.8660254037844386im
            @test greens_function_bethe_analytic(0, 10) == -0.2im
            # Vector{Complex}
            g = greens_function_bethe_analytic([2.0 + 0.1im, 3.0 + 0.5im])
            @test typeof(g) === Vector{ComplexF64}
            @test length(g) === 2
            @test g[1] == greens_function_bethe_analytic(2.0 + 0.1im)
            @test g[2] == greens_function_bethe_analytic(3.0 + 0.5im)
        end # analytic

        @testset "simple" begin
            # 101 poles
            G = greens_function_bethe_simple(101)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test length(G) === 101
            @test moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test location(G, 51) ≈ 0 atol = 10 * eps()
            @test norm(locations(G) + reverse(locations(G))) < 50 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 600 * eps()
            # 100 poles
            G = greens_function_bethe_simple(100)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test length(G) === 100
            @test moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test norm(locations(G) + reverse(locations(G))) < 100 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 600 * eps()
            # 101 poles, D = 2
            G = greens_function_bethe_simple(101, 2)
            @test moment(G, 0) ≈ 1 rtol = 10 * eps()
        end # simple

        @testset "grid" begin
            # 1 pole
            W = [0.0]
            G = greens_function_bethe_grid(W)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test isone(length(G))
            @test locations(G) == W
            @test locations(G) !== W
            @test only(weights(G)) === 1.0
            # 101 poles
            W = range(-1, 1; length = 101)
            G = greens_function_bethe_grid(W)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test length(G) === 101
            @test locations(G) == W
            @test locations(G) !== W
            @test moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 10 * eps()
            @test weight(G, 51) ≈ 0.012732183237577577 atol = eps()
            # 100 poles
            W = range(-1, 1; length = 100)
            G = greens_function_bethe_grid(W)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test length(G) === 100
            @test locations(G) == W
            @test locations(G) !== W
            @test moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 10 * eps()
            @test weight(G, 51) ≈ 0.012860130639746004 atol = eps()
            # 101 poles, D = 2
            W = range(-3, 3; length = 101)
            G = greens_function_bethe_grid(W, 2)
            @test moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 10 * eps()
            @test all(iszero, view(weights(G), 1:17))
            @test all(iszero, view(weights(G), 85:101))
            @test weight(G, 51) ≈ 0.01909787694960996 atol = eps()
            # non-equidistant grid
            # Test if dense grid in middle has smaller weights.
            W = [-1:0.01:-0.51; -0.5:0.005:0.5; 0.51:0.01:1]
            G = greens_function_bethe_grid(W)
            w1 = weight(G, 50)
            @test all(i -> i < w1, view(weights(G), 51:251))
            @test weight(G, 151) ≈ 0.0031830955461067956 atol = eps()
        end # grid

        @testset "grid Hubbard III" begin
            # 1 pole
            G = greens_function_bethe_grid_hubbard3([5.0])
            @test locations(G) == [5.0]
            @test weights(G) == [1.0]
            # uniform grid
            grid = range(-5, 5; length = 101)
            # U = 0
            G = greens_function_bethe_grid_hubbard3(grid)
            G0 = greens_function_bethe_grid(grid)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test length(G) === 101
            @test locations(G) == grid
            @test locations(G) !== grid
            @test norm(weights(G) - weights(G0)) < 10 * eps()
            # U = 3
            G = greens_function_bethe_grid_hubbard3(grid, 3)
            @test amplitudes(G)[36] ≈ 0.1783752245364157 atol = 10 * eps()
            @test amplitudes(G)[51] == 0
            @test amplitudes(G)[66] ≈ 0.1783752245364157 atol = 10 * eps()
            @test moment(G, 0) ≈ 1 atol = 10 * eps()
        end # grid Hubbdard III

        @testset "equal weight" begin
            @test_throws DomainError greens_function_bethe_equal_weight(2)
            G = greens_function_bethe_equal_weight(101)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test length(G) === 101
            @test all(i -> i === 1 / 101, weights(G))
            @test norm(locations(G) + reverse(locations(G))) === 0.0
        end # equal weight
    end # Bethe lattice

    @testset "user supplied dispersion" begin
        H_k = [[1 + 0im 2; 2 1], [3 4; 4 3]]
        Σ_stat = Diagonal([1, 0])
        Σ_dyn = PolesSumBlock([-2, 3], [[0 im; -im 5], [0 0; 0 6]]) # self-energy only on [2, 2] index

        @testset "non-interacting" begin
            @inferred greens_function_local(H_k, 0.5)
            G0 = greens_function_local(H_k, 0.5)
            @test length(G0) == 3
            @test locations(G0) == [-1, 3, 7] .- 0.5
            @test norm(moment(G0, 0) - I) < 10 * eps() # normalized
            @test norm(weight(G0, 1) - [0.5 -0.5; -0.5 0.5]) < 2 * eps()
            @test norm(weight(G0, 2) - [0.25 0.25; 0.25 0.25]) < 2 * eps()
            @test norm(weight(G0, 3) - [0.25 0.25; 0.25 0.25]) < 2 * eps()
        end # non-interacting

        @testset "interacting" begin
            # μ = 0
            @inferred greens_function_local(H_k, Σ_stat, Σ_dyn, 0)
            G = greens_function_local(H_k, Σ_stat, Σ_dyn, 0)
            locs_ref = [
                -3.672298697206984,
                -3.4963335910609388,
                -0.13601938094586696,
                -0.11866710984063822,
                2.4215987184100305,
                3.2883126139747096,
                5.386719359742839,
                8.326688086926891,
            ]
            @test all(<(10 * eps()), locations(G) - locs_ref)
            @test norm(moment(G, 0) - I) < 10 * eps()
            m1 = [3 3; 3 2]
            @test norm(moment(G, 1) - m1) < 100 * eps()

            # μ = 4 must not be a simple shift of poles
            G = greens_function_local(H_k, Σ_stat, Σ_dyn, 4)
            locs_ref = [
                -5.966810528426411,
                -5.796191057467859,
                -2.0437833375992263,
                -1.5551427327462921,
                -0.06386413066063401,
                2.0401530820076825,
                4.074457996686286,
                5.311180708206492,
            ]
            @test all(<(10 * eps()), locations(G) - locs_ref)
            @test norm(moment(G, 0) - I) < 10 * eps()
            m1 = [-1 3; 3 -2]
            @test norm(moment(G, 1) - m1) < 100 * eps()
        end # interacting
    end # user supplied dispersion
end # Green's function
