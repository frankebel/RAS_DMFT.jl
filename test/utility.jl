using Fermions
using Fermions.Lanczos
using Fermions.Wavefunctions
using LinearAlgebra
using RAS_DMFT
using Random
using Test

@testset "util" begin
    @testset "get_CI_parameters" begin
        @test get_CI_parameters(10, 5, 1, 1) == (4, 3, 3)
        @test get_CI_parameters(10, 6, 1, 1) == (4, 4, 2)
        @test get_CI_parameters(10, 4, 1, 1) == (4, 2, 4)
        @test get_CI_parameters(10, 5, 2, 1) == (5, 2, 3)
        @test get_CI_parameters(10, 5, 1, 2) == (5, 3, 2)
        @test get_CI_parameters(10, 6, 2, 1) == (5, 3, 2)
        @test get_CI_parameters(10, 4, 1, 2) == (5, 2, 3)
        # non-sensical input values still work
        @test get_CI_parameters(10, 10, 1, 1) == (4, 8, -2)
        @test get_CI_parameters(10, 0, 1, 1) == (4, -2, 8)
    end # get_CI_parameters

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

    @testset "δ_gaussian" begin
        @test δ_gaussian(0.01, 0.04, 1.0, 0.0) ≈ 0.01 atol = 2 * eps() # ω = 0
        @test δ_gaussian(0.01, 0.04, 1.0, 10.0) ≈ 0.04 atol = 2 * eps() # ω → ∞
        @test δ_gaussian(0.01, 0.04, 1.0, 0.5) === δ_gaussian(0.01, 0.04, 1.0, -0.5) # symmetry
        @test δ_gaussian(0.01, 0.04, 1.0, 0.2) < δ_gaussian(0.01, 0.04, 1.0, 0.3)
    end # δ_gaussian

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

        # finite broadening
        # Create a symmetric uniform density on `n_tot` poles.
        n_tot = 100 # number of poles ≙ filling for μ = ∞
        h = Diagonal(range(-1, 1, n_tot)) # uniform density in [-1, 1]
        Hk = [h]
        Z = (-10:0.01:10) .+ 0.05im
        Σ = [zero(h) for _ in eachindex(Z)]
        # test half-filling
        μ, filling = find_chemical_potential(Z, Hk, Σ, n_tot / 2)
        @test μ ≈ 0 atol = 2.0e-3
        @test filling ≈ n_tot / 2 atol = 2.0e-2
    end # find chemical potential

    @testset "_issorted_and_unique" begin
        @test RAS_DMFT._issorted_and_unique(1:10)
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([1, 0]) # not sorted
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([0, 0, 1]) # not unique
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([-0.0, 0.0]) # duplicate zeros
    end # _issorted_and_unique

    @testset "moment" begin
        ref = [0.1, 0.2, 0.3, 0.4, 0.5] # must not matter
        imf = [-0.2, -0.4, 0.0, -0.6, -0.8] # asymmetric input
        f = ref + im * imf
        W = [-0.4, -0.1, 0.0, 0.5, 1.0] # non-equidistant grid
        @test moment(f, W) == 0.61 / π
        @test moment(f, W, 0) == 0.61 / π
        @test moment(f, W, 1) == 0.33 / π
        @test moment(f, W, 2) == 0.2806 / π
    end # moment

    @testset "moments" begin
        ref = [0.1, 0.2, 0.3, 0.4, 0.5] # must not matter
        imf = [-0.2, -0.4, 0.0, -0.6, -0.8] # asymmetric input
        f = ref + im * imf
        W = [-0.4, -0.1, 0.0, 0.5, 1.0] # non-equidistant grid
        @test moments(f, W, 0:2) == [0.61, 0.33, 0.2806] ./ π
    end # moments

    @testset "quasiparticle_weight" begin
        # pole
        Σ = PolesSum([-0.25, -0.01, 0.5], [1.0, 2.0, 3.0])
        @test quasiparticle_weight(Σ) == inv(20029)
        @test quasiparticle_weight(Σ; tol = 1.0) == inv(20029)
        @test quasiparticle_weight(Σ; tol = 1.1) == inv(20013)
        @test quasiparticle_weight(Σ; λ = 1.0e-2) ≈ inv(10028.9696428138) atol = 10 * eps()
        @test quasiparticle_weight(Σ; tol = 1.1, λ = 1.0e-2) ≈ inv(10012.995201919231) atol = 10 * eps()

        # grid
        W = [-0.5, -0.25, 0.0, 0.25, 0.5]
        Σ = [2.0 + 0.1im, 1.5 + 0.2im, 1.0 + 0.3im, 0.5 + 0.4im, 0.0 + 0.5im]
        @test quasiparticle_weight(Σ, W, 1) == inv(1 + 2)
    end # quasiparticle_weight

    @testset "quasiparticle_weight_gaussian" begin
        Σ = PolesSum([-0.25, -0.01, 0.5], [1.0, 2.0, 3.0])
        @test quasiparticle_weight_gaussian(Σ, 0.01, 0.01) ≈ -0.00015696719254258983 atol =
            100 * eps()
        @test quasiparticle_weight_gaussian(Σ, 0.001, 0.01) ≈ -0.00018214101681205788 atol =
            100 * eps()
    end # quasiparticle_weight_gaussian
end # util
