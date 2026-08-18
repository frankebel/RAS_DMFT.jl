using Fermions
using Fermions.Lanczos
using Fermions.Wavefunctions
using LinearAlgebra
using RAS_DMFT
using Random
using Test

@testset "util" begin
    @testset "get_RAS_parameters" begin
        @test get_RAS_parameters(10, 5, 1, 1) == (4, 3, 3)
        @test get_RAS_parameters(10, 6, 1, 1) == (4, 4, 2)
        @test get_RAS_parameters(10, 4, 1, 1) == (4, 2, 4)
        @test get_RAS_parameters(10, 5, 2, 1) == (5, 2, 3)
        @test get_RAS_parameters(10, 5, 1, 2) == (5, 3, 2)
        @test get_RAS_parameters(10, 6, 2, 1) == (5, 3, 2)
        @test get_RAS_parameters(10, 4, 1, 2) == (5, 2, 3)
        # non-sensical input values still work
        @test get_RAS_parameters(10, 10, 1, 1) == (4, 8, -2)
        @test get_RAS_parameters(10, 0, 1, 1) == (4, -2, 8)
    end # get_RAS_parameters

    @testset "init_system" begin
        # parameters
        n_bath = 31
        U = 4.0
        μ = U / 2
        L_v = 1
        L_c = 1
        p = 2
        var = eps()

        E0_target = -21.527949990417255 # target ground state energy
        Δ = hybridization_function_bethe_simple(n_bath)
        fs = FockSpace(Orbitals(2 + L_v + L_c), FermionicSpin(1 // 2))
        n = occupations(fs)
        H_int = U * n[1, -1 // 2] * n[1, 1 // 2]
        H, E0, ψ0 = init_system(Δ, H_int, -μ, L_v, L_c, p, eps())
        Hψ = H * ψ0
        variance = Hψ ⋅ Hψ
        @test variance < var
        @test E0 ≈ E0_target rtol = 2.0e-13
    end # init system

    @testset "Kondo temperature" begin
        @test temperature_kondo(0.3, -0.1, 0.1) == 0.04297872341114842
        @test temperature_kondo(0.2, -0.1, 0.015) == 0.00020610334475146955
    end # Kondo temperature

    @testset "find chemical potential" begin
        Random.seed!(0)
        n_b = 4
        n_fill = n_b * 0.7
        H_k = [Hermitian(randn(n_b, n_b)) for _ in 1:10]
        Σ_stat = Diagonal(randn(n_b))
        amps = [randn(4, 2) for _ in 1:20]
        wgts = [Hermitian(amp * amp') for amp in amps]
        Σ_dyn = PolesSumBlock(randn(20) .* 2, wgts)
        μ, n = find_chemical_potential(H_k, Σ_stat, Σ_dyn, n_fill; μ_min = -1, μ_max = 8)
        @test μ ≈ 5.823782742023468 atol = 1.0e-6
        @test n ≈ n_fill atol = 1.0e-7
    end # find chemical potential

    @testset "_issorted_and_unique" begin
        @test RAS_DMFT._issorted_and_unique(1:10)
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([1, 0]) # not sorted
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([0, 0, 1]) # not unique
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([-0.0, 0.0]) # duplicate zeros
    end # _issorted_and_unique

    @testset "quasiparticle_weight" begin
        Σ = PolesSum([-0.25, -0.01, 0.5], [1.0, 2.0, 3.0])
        @test quasiparticle_weight(Σ) == inv(20029)
        @test quasiparticle_weight(Σ; tol = 1.0) == inv(20029)
        @test quasiparticle_weight(Σ; tol = 1.1) == inv(20013)
        @test quasiparticle_weight(Σ; λ = 1.0e-2) ≈ inv(10028.9696428138) atol = 10 * eps()
        @test quasiparticle_weight(Σ; tol = 1.1, λ = 1.0e-2) ≈ inv(10012.995201919231) atol = 10 * eps()
    end # quasiparticle_weight

    @testset "quasiparticle_weight_inflections" begin
        # single pole with optimum λ = 1
        Σ = PolesSum([1.0], [2.0])
        @test quasiparticle_weight_inflections(Σ; λmax = 2.0) ≈ [1.0] atol = 2 * eps()
        # outside the search window
        @test isempty(quasiparticle_weight_inflections(Σ; λmax = 0.5))
        @test isempty(quasiparticle_weight_inflections(Σ; λmin = 1.5, λmax = 2.0))
        @test isempty(quasiparticle_weight_inflections(Σ; tol = 3.0))
        # invalid arguments
        @test_throws ArgumentError quasiparticle_weight_inflections(Σ; tol = -1)
        @test_throws ArgumentError quasiparticle_weight_inflections(Σ; λmin = 0)
        @test_throws ArgumentError quasiparticle_weight_inflections(Σ; λmax = 0)

        Σ2 = PolesSum([-1.0, 1.0], [0.5, 0.5])
        @test quasiparticle_weight_inflections(Σ2; λmax = 2.0) ≈ [sqrt(2 / 3)] atol = 2 * eps()

        #  three inflections
        Σ3 = PolesSum([1.0e-4, 1.0], [1.0e-7, 2.0])
        λ3 = quasiparticle_weight_inflections(Σ3; λmax = 2.0)
        @test λ3 ≈ [1.2018504652954801e-4, 1.9681873959735308e-2, 0.9999998333332845] atol = 10 * eps()
    end # quasiparticle_weight_inflections
end # util
