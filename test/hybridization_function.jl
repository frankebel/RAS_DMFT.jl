using RAS_DMFT
using LinearAlgebra
using Test

@testset "hybridization" begin
    @testset "Bethe lattice" begin
        @testset "analytic" begin
            # test different number types
            @test hybridization_function_bethe_analytic(-2) == -0.1339745962155614
            @test hybridization_function_bethe_analytic(-true) == -1 // 2
            @test hybridization_function_bethe_analytic(-1 // 2) ==
                -0.25 - 0.4330127018922193im
            @test hybridization_function_bethe_analytic(false) == -0.5im
            @test hybridization_function_bethe_analytic(0.1im) == -0.4524937810560445im
            @test hybridization_function_bethe_analytic(0.5) == 0.25 - 0.4330127018922193im
            @test hybridization_function_bethe_analytic(0x01) == 0.5
            @test hybridization_function_bethe_analytic(2.0) == 0.1339745962155614
            # # vary half-bandwidth
            @test hybridization_function_bethe_analytic(-1, 2) ==
                -0.5 - 0.8660254037844386im
            @test hybridization_function_bethe_analytic(0, 10) == -5im
            # Vector{Complex}
            Δ = hybridization_function_bethe_analytic([2.0 + 0.1im, 3.0 + 0.5im])
            @test typeof(Δ) === Vector{ComplexF64}
            @test length(Δ) === 2
            @test Δ[1] == hybridization_function_bethe_analytic(2.0 + 0.1im)
            @test Δ[2] == hybridization_function_bethe_analytic(3.0 + 0.5im)
        end # analytic

        @testset "simple" begin
            # 101 poles
            Δ = hybridization_function_bethe_simple(101)
            @test typeof(Δ) == PolesSum{Float64, Float64}
            @test length(Δ) === 101
            @test moment(Δ, 0) ≈ 0.25 rtol = 10 * eps()
            @test location(Δ, 51) ≈ 0 atol = 10 * eps()
            @test norm(locations(Δ) + reverse(locations(Δ))) < 50 * eps()
            @test norm(weights(Δ) - reverse(weights(Δ))) < 600 * eps()
            # 100 poles
            Δ = hybridization_function_bethe_simple(100)
            @test typeof(Δ) === PolesSum{Float64, Float64}
            @test length(Δ) === 100
            @test moment(Δ, 0) ≈ 0.25 rtol = 10 * eps()
            @test norm(locations(Δ) + reverse(locations(Δ))) < 100 * eps()
            @test norm(weights(Δ) - reverse(weights(Δ))) < 600 * eps()
            # 101 poles, D = 2
            Δ = greens_function_bethe_simple(101, 2)
            @test moment(Δ, 0) ≈ 1.0 rtol = 10 * eps()
        end # simple

        @testset "grid Hubbard III" begin
            grid = range(-5, 5; length = 101)
            # U = 0
            Δ = hybridization_function_bethe_grid_hubbard3(grid)
            Δ0 = hybridization_function_bethe_grid(grid)
            @test typeof(Δ) === PolesSum{Float64, Float64}
            @test length(Δ) === 101
            @test locations(Δ) == grid
            @test locations(Δ) !== grid
            @test norm(amplitudes(Δ) - amplitudes(Δ0)) < 10 * eps()
            # U = 3
            Δ = hybridization_function_bethe_grid_hubbard3(grid, 3)
            @test amplitudes(Δ)[36] ≈ 0.08918761226820784 atol = 10 * eps()
            @test amplitudes(Δ)[51] == 0
            @test amplitudes(Δ)[66] ≈ 0.08918761226820784 atol = 10 * eps()
            @test moment(Δ, 0) ≈ 0.25 atol = 10 * eps()
        end # grid Hubbard III

        @testset "grid" begin
            # 101 poles
            W = range(-1, 1; length = 101)
            Δ = hybridization_function_bethe_grid(W)
            @test typeof(Δ) === PolesSum{Float64, Float64}
            @test length(Δ) === 101
            @test locations(Δ) == W
            @test locations(Δ) !== W
            @test moment(Δ, 0) ≈ 0.25 rtol = 10 * eps()
            @test norm(weights(Δ) - reverse(weights(Δ))) < 10 * eps()
            @test amplitudes(Δ)[51] ≈ 0.056418488187777546 atol = eps()
            # 100 poles
            W = range(-1, 1; length = 100)
            Δ = hybridization_function_bethe_grid(W)
            @test typeof(Δ) === PolesSum{Float64, Float64}
            @test length(Δ) === 100
            @test locations(Δ) == W
            @test locations(Δ) !== W
            @test moment(Δ, 0) ≈ 0.25 rtol = 10 * eps()
            @test norm(weights(Δ) - reverse(weights(Δ))) < 10 * eps()
            @test amplitudes(Δ)[51] ≈ 0.05670125801017559 atol = eps()
            # 101 poles, D = 2
            W = range(-3, 3; length = 101)
            Δ = hybridization_function_bethe_grid(W, 2)
            @test sum(weights(Δ)) ≈ 1 rtol = 10 * eps()
            @test norm(weights(Δ) - reverse(weights(Δ))) < 10 * eps()
            @test all(iszero, view(amplitudes(Δ), 1:17))
            @test all(iszero, view(amplitudes(Δ), 85:101))
            @test amplitudes(Δ)[51] ≈ 0.13819506847065838 atol = eps()
            # non-equidistant grid
            # Test if dense grid in middle has smaller weights.
            W = [-1:0.01:-0.51; -0.5:0.005:0.5; 0.51:0.01:1]
            Δ = hybridization_function_bethe_grid(W)
            w1 = amplitudes(Δ)[50]
            @test all(i -> i < w1, view(amplitudes(Δ), 51:251))
            @test amplitudes(Δ)[151] ≈ 0.028209464484933045 atol = eps()
        end # grid

        @testset "equal weight" begin
            @test_throws DomainError hybridization_function_bethe_equal_weight(2)
            # D = 1
            Δ = hybridization_function_bethe_equal_weight(101)
            @test typeof(Δ) === PolesSum{Float64, Float64}
            @test length(Δ) === 101
            @test all(i -> i === 1 / sqrt(101) / 2, amplitudes(Δ))
            @test norm(locations(Δ) + reverse(locations(Δ))) === 0.0
            # D = 4
            Δ = hybridization_function_bethe_equal_weight(101, 4)
            @test all(i -> i === 2 / sqrt(101), amplitudes(Δ))
        end # equal weight
    end # Bethe lattice
end # hybridization
