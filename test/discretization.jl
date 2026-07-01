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
    H, E0, ψ0 = init_system(Δ0, H_int, ϵ_imp, n_v_bit, n_c_bit, e, var)

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

        # Lorentzian
        # self-energy
        Σ = self_energy_IFG_lorentzian(Σ_H, C, W, δ)
        # new hybridization function
        Δ_grid = Δ0_analytic(Z .+ μ - Σ)
        Δ_new_L = equal_weight_discretization(-imag(Δ_grid), W, δ, n_bath)
        @test typeof(Δ_new_L) === PolesSum{Float64, Float64}
        @test typeof(Δ_grid) === Vector{ComplexF64}
        @test length(Δ_new_L) === n_bath
        @test all(b -> isapprox(b, 1 / sqrt(n_bath) / 2; rtol = 3.0e-3), amplitudes(Δ_new_L))
        # small weight loss due to truncated interval
        @test 0.2486 <= moment(Δ_new_L, 0) <= 0.25
        @test moment(Δ_new_L, 1) < 200 * eps() # PHS

        # Gaussian
        # self-energy
        Σ = self_energy_IFG_gaussian(Σ_H, C, W, δ)
        # new hybridization function
        Δ_grid = Δ0_analytic(Z .+ μ - Σ)
        Δ_new_G = equal_weight_discretization(-imag(Δ_grid), real(Z), δ, n_bath)
        @test typeof(Δ_new_G) === PolesSum{Float64, Float64}
        @test typeof(Δ_grid) === Vector{ComplexF64}
        @test length(Δ_new_G) === n_bath
        weights_without_zero = copy(weights(Δ_new_G))
        popat!(weights_without_zero, cld(n_bath, 2))
        @test all(b -> isapprox(b, inv(4 * n_bath); atol = 2.0e-5), weights_without_zero)
        @test weight(Δ_new_G, cld(n_bath, 2)) ≈ inv(4 * n_bath) atol = 2.0e-4
        # small weight loss due to truncated interval
        @test 0.2486 <= moment(Δ_new_G, 0) <= 0.25
        @test moment(Δ_new_G, 1) < 200 * eps() # PHS

        # must not be equal
        @test locations(Δ_new_L) != locations(Δ_new_G)
        @test amplitudes(Δ_new_L) != amplitudes(Δ_new_G)
    end # block Lanczos

    # Discretization for Gaussian returned wrong number of poles.
    w = -10:0.0002:10
    g = similar(w)
    @. g = exp(-w^2)
    wrong_length = 0
    for i in 3:2:1001
        # test different number of poles
        dis = equal_weight_discretization(g, w, 0.04, i)
        length(dis) == i || (wrong_length += 1)
    end
    @test iszero(wrong_length)

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
