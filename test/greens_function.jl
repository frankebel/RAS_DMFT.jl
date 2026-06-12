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
            @test RAS_DMFT.moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test locations(G)[51] ≈ 0 atol = 10 * eps()
            @test norm(locations(G) + reverse(locations(G))) < 50 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 600 * eps()
            # 100 poles
            G = greens_function_bethe_simple(100)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test length(G) === 100
            @test RAS_DMFT.moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test norm(locations(G) + reverse(locations(G))) < 100 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 600 * eps()
            # 101 poles, D = 2
            G = greens_function_bethe_simple(101, 2)
            @test RAS_DMFT.moment(G, 0) ≈ 1 rtol = 10 * eps()
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
            @test RAS_DMFT.moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 10 * eps()
            @test weight(G, 51) ≈ 0.012732183237577577 atol = eps()
            # 100 poles
            W = range(-1, 1; length = 100)
            G = greens_function_bethe_grid(W)
            @test typeof(G) === PolesSum{Float64, Float64}
            @test length(G) === 100
            @test locations(G) == W
            @test locations(G) !== W
            @test RAS_DMFT.moment(G, 0) ≈ 1 rtol = 10 * eps()
            @test norm(weights(G) - reverse(weights(G))) < 10 * eps()
            @test weight(G, 51) ≈ 0.012860130639746004 atol = eps()
            # 101 poles, D = 2
            W = range(-3, 3; length = 101)
            G = greens_function_bethe_grid(W, 2)
            @test RAS_DMFT.moment(G, 0) ≈ 1 rtol = 10 * eps()
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
            @test RAS_DMFT.moment(G, 0) ≈ 1 atol = 10 * eps()
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
        Hk = [[1 + 0.0im 2; 2 1], [3 4; 4 3]]
        Z = (-10:-9) .+ 0.1im
        Σ = [Diagonal([0, 5 + im]), Diagonal([0, 6 + im])] # self-energy only on [2, 2] index

        @testset "non-interacting" begin
            @inferred greens_function_local(Hk)
            G0 = greens_function_local(Hk)
            @test length(G0) == 3
            @test locations(G0) == [-1, 3, 7]
            @test norm(RAS_DMFT.moment(G0, 0) - I) < 10 * eps() # normalized
            @test norm(weight(G0, 1) - [0.5 -0.5; -0.5 0.5]) < 2 * eps()
            @test norm(weight(G0, 2) - [0.25 0.25; 0.25 0.25]) < 2 * eps()
            @test norm(weight(G0, 3) - [0.25 0.25; 0.25 0.25]) < 2 * eps()

            # finite broadening
            G0 = greens_function_local(Z, 0, Hk)
            @test length(G0) == 2
            @test norm(
                G0[1] - [
                    -0.08948370259088873 - 0.0008516301906910084im 0.021613692792397263 + 0.0003827853135677249im
                    0.021613692792397263 + 0.0003827853135677249im -0.08948370259088873 - 0.0008516301906910084im
                ],
            ) < eps()
            @test norm(
                G0[2] - [
                    -0.09894651224745543 - 0.001052379439830884im 0.0260339595538256 + 0.0005098764576851289im
                    0.0260339595538256 + 0.0005098764576851289im -0.09894651224745543 - 0.001052379439830884im
                ],
            ) < eps()
        end # non-interacting

        @testset "interacting" begin
            G = greens_function_local(Z, 0, Hk, Σ)
            @test length(G) == 2
            @test norm(
                G[1] - [
                    -0.08778121899483271 - 0.0005616185151491649im 0.014949094763816347 - 0.0006950450017352803im
                    0.014949094763816347 - 0.0006950450017352803im -0.061604402244826835 + 0.0034066891091674998im
                ],
            ) < eps()
            @test norm(
                G[2] - [
                    -0.09626381282560724 - 0.0006774835502274772im 0.01636751358241116 - 0.0007517320959592241im
                    0.01636751358241116 - 0.0007517320959592241im -0.06186055815770346 + 0.003430278158287015im
                ],
            ) < eps()
        end # interacting

        @testset "spectrum Gauss" begin
            W = [-1, 1]
            A = spectral_function_gauss(W, 0, Hk, 0.05)
            @test typeof(A) === Vector{Matrix{Float64}}
            @test norm(
                A[1] - [
                    3.989422804014326 -3.989422804014326
                    -3.989422804014326 3.989422804014326
                ],
            ) < eps()
            @test iszero(A[2])
        end # spectrum Gauss
    end # user supplied dispersion
end # Green's function
