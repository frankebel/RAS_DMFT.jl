using RAS_DMFT
using Fermions
using LinearAlgebra
using Test

@testset "self-energy" begin
    # parameters
    n_bath = 31
    U = 4.0
    μ = U / 2
    ϵ_imp = -μ
    L_v = 1
    L_c = 1
    p = 2
    n_kryl = 100
    var = eps()
    tol = 1.0e-8 # merge weights smaller than this

    Δ0 = hybridization_function_bethe_simple(n_bath)
    # Operators for positive frequencies. Negative ones are calculated by adjoint.
    fs = FockSpace(Orbitals(2 + L_v + L_c), FermionicSpin(1 // 2))
    c = annihilators(fs)
    n = occupations(fs)
    d_dag = c[1, -1 // 2]' # d_↓^†

    # only interacting part
    H_int = U * n[1, 1 // 2] * n[1, -1 // 2]
    q_dag = H_int * d_dag - d_dag * H_int  # q_↓^† = [H_int, d^†]
    H, _, ψ0 = init_system(Δ0, H_int, ϵ_imp, 0, L_v, L_c, p, var)
    O_Σ_H = q_dag' * d_dag + d_dag * q_dag'
    Σ_H = dot(ψ0, O_Σ_H, ψ0)

    # linear shift
    q_dag_tilde = q_dag - Σ_H * d_dag
    O = [q_dag_tilde, d_dag]

    # impurity solver
    C_plus = correlator_plus(H, ψ0, O, n_kryl)
    C_minus = correlator_minus(H, ψ0, map(adjoint, O), n_kryl)
    C = transpose(C_minus) + C_plus
    merge_small_weight!(C, tol)

    @testset "Dyson" begin
        # impurity Green's function
        G_plus = PolesSum(C_plus, 2, 2)
        remove_zero_weight!(G_plus)
        merge_degenerate_poles!(G_plus)
        merge_small_weight!(G_plus, 1.0e-11)
        G_minus = flip_spectrum(G_plus)
        G_imp = G_minus + G_plus

        Σ_H, Σ = self_energy_dyson(-μ, Δ0, G_imp, -5:0.02:5)
        @test Σ_H ≈ U / 2 atol = 10 * eps() # half-filling
        @test moment(Σ, 0) ≈ U^2 / 4 atol = 2.0e-11
        @test !any(iszero, locations(Σ)) # no pole at 0 for metal
    end # self_energy_dyson

    @testset "Schur" begin
        Σ = PolesSum(self_energy_schur(C), 1, 1)
        merge_small_weight!(Σ, tol)
        @test moment(Σ, 0) ≈ U^2 / 4 rtol = 80 * eps()
        @test moment(Σ, 1) ≈ 0 atol = 5.0e-10
    end # Schur

    @testset "Schur block size four" begin
        locsC = [-1.5, -0.5, 0.8]
        # lower-right 2x2 blocks are dyadic and sum is exactly I
        W1 = [1.2 0.3 0.1 0.0; 0.3 0.9 0.2 0.1; 0.1 0.2 0.25 0.125; 0.0 0.1 0.125 0.25]
        W2 = [0.8 0.2 0.0 0.2; 0.2 1.0 0.1 0.0; 0.0 0.1 0.375 0.0; 0.2 0.0 0.0 0.375]
        W3 = [1.4 0.1 0.2 0.3; 0.1 0.7 0.0 0.1; 0.2 0.0 0.375 -0.125; 0.3 0.1 -0.125 0.375]
        C4 = PolesSumBlock(locsC, [W1, W2, W3])

        # assert desired properties
        m0 = moment(C4, 0)
        @test m0[3:4, 3:4] == I # Green's function moment from anti-commutator
        @test norm(m0[1:2, 3:4]) > 0.1 * norm(m0) # coupling must not vanish

        Σ = self_energy_schur(C4)
        merge_degenerate_poles!(Σ, 30 * eps())
        @test size(Σ) == (2, 2)
        @test length(Σ) == 15
        # zeroth moment: invert, project, invert
        @test moment(Σ, 0) ≈ inv((inv(m0))[1:2, 1:2]) atol = 1.0e2 * eps()
        # m1, m2 not analytically tested
        @test moment(Σ, 1) ≈ [-1.4015 -0.47125; -0.47125 -1.120625] atol = 1.0e3 * eps()
        @test moment(Σ, 2) ≈
            [3.5254375 0.58536875; 0.58536875 2.412646875] atol = 1.0e3 * eps()
    end # Schur block size four
end # self-energy
