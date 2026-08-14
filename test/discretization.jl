using RAS_DMFT
using Fermions
using LinearAlgebra
using Test

@testset "discretization" begin
    # parameters
    n_bath = 101
    U = 2.0
    μ = U / 2
    ϵ_imp = -μ
    n_v_bit = 1
    n_c_bit = 1
    e = 1
    n_kryl = 100
    var = eps()
    W = -5:0.001:5
    δ = 0.04
    # do not change parameters below
    Z = W .+ im * δ

    Δ0_analytic = hybridization_function_bethe_analytic
    Δ0 = hybridization_function_bethe_simple(n_bath)

    # Operators for positive frequencies. Negative ones are calculated by adjoint.
    fs = FockSpace(Orbitals(2 + n_v_bit + n_c_bit), FermionicSpin(1 // 2))
    c = annihilators(fs)
    n = occupations(fs)
    H_int = U * n[1, 1 // 2] * n[1, -1 // 2]
    d_dag = c[1, -1 // 2]' # d_↓^†
    q_dag = H_int * d_dag - d_dag * H_int  # q_↓^† = [H_int, d^†]
    O = [q_dag, d_dag]

    # initialize system
    H, _, ψ0 = init_system(Δ0, H_int, ϵ_imp, n_v_bit, n_c_bit, e, var)

    @testset "Lanczos" begin
        # impurity Green's functions
        G_plus = correlator_plus(H, ψ0, d_dag, n_kryl)
        G_minus = correlator_minus(H, ψ0, d_dag', n_kryl)
        G_imp = G_plus + G_minus

        # self-energy
        Σ_H, Σ = self_energy_dyson(ϵ_imp, Δ0, G_imp, W)
        merge_small_weight!(Σ, 1.0e-11)

        # new hypridization function
        Δ = update_hybridization_function(Δ0, μ, Σ_H, Σ)
        merge_degenerate_poles!(Δ)
        merge_small_weight!(Δ, 1.0e-11)
        Δ_new = discretize_similar_weight(Δ, sqrt(eps()), n_bath)

        @test length(Δ_new) == 101
        @test iszero(location(Δ_new, 51))
        @test weight(Δ_new, 51) ≈ 0.0020833288949242113 atol = 1.0e-7
        @test moment(Δ_new, 0) ≈ 0.25 atol = 10 * eps()
        @test moment(Δ_new, 1) ≈ 0.0 atol = sqrt(eps())

        # insulating solution
        P = PolesSum([rand(100) .- 2; rand(100) .+ 1], rand(200))
        sort!(P)
        foo = discretize_similar_weight(P, 0.01, 11)
        @test iszero(location(foo, 6))
        @test iszero(weight(foo, 6))
    end # Lanczos

    @testset "block Lanczos" begin
        # Hartree term
        O_H = O[1]' * O[2] + O[2] * O[1]'
        Σ_H = dot(ψ0, O_H, ψ0)

        # impurity correlators
        C_plus = correlator_plus(H, ψ0, O, n_kryl)
        C_minus = correlator_minus(H, ψ0, map(adjoint, O), n_kryl)
        C = transpose(C_minus) + C_plus
    end # block Lanczos

    @testset "discretize_to_grid" begin
        W = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        f = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0] # not same step size
        # 1 pole
        grid = [0.0]
        foo = discretize_to_grid(f, W, grid)
        @test locations(foo) == grid
        @test locations(foo) !== grid
        @test weights(foo) == [15.125 / π]
        # 1 pole < W
        foo = discretize_to_grid(f, W, [-10.0])
        @test locations(foo) == [-10.0]
        @test weights(foo) == [15.125 / π]
        # 1 pole > W
        foo = discretize_to_grid(f, W, [10.0])
        @test locations(foo) == [10.0]
        @test weights(foo) == [15.125 / π]
        # 2 asymmetric poles
        foo = discretize_to_grid(f, W, [-3.0, 5.0])
        @test locations(foo) == [-3.0, 5.0]
        @test weights(foo) == [6.0 / π, 9.125 / π]
    end # discretize_to_grid
end # discretization
