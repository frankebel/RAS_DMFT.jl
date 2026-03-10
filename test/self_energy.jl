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
    step_size = 0.02
    W = -10:step_size:10
    δ = 0.08
    tol = 1.0e-8 # merge weights smaller than this

    Δ0 = hybridization_function_bethe_simple(n_bath)
    # Operators for positive frequencies. Negative ones are calculated by adjoint.
    fs = FockSpace(Orbitals(2 + L_v + L_c), FermionicSpin(1 // 2))
    c = annihilators(fs)
    n = occupations(fs)
    H_int = U * n[1, 1 // 2] * n[1, -1 // 2]
    d_dag = c[1, -1 // 2]' # d_↓^†

    # only interacting part
    H_int = U * n[1, 1 // 2] * n[1, -1 // 2]
    q_dag = H_int * d_dag - d_dag * H_int  # q_↓^† = [H_int, d^†]
    H, E0, ψ0 = init_system(Δ0, H_int, ϵ_imp, L_v, L_c, p, var)
    O_Σ_H = q_dag' * d_dag + d_dag * q_dag'
    Σ_H = dot(ψ0, O_Σ_H, ψ0)

    # quadratic shift
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
        @test Σ_H ≈ U / 2 atol = 100 * eps() # half-filling
        @test RAS_DMFT.moment(Σ, 0) ≈ U^2 / 4 atol = 1.0e-5 # bad agreement
        @test !any(iszero, locations(Σ)) # no pole at 0 for metal
    end # self_energy_dyson

    @testset "IFG" begin
        # on the real axis
        Σ = PolesSum(self_energy_IFG(C), 1, 1)
        merge_small_weight!(Σ, tol)
        @test RAS_DMFT.moment(Σ, 0) ≈ U^2 / 4 rtol = 1.0e3 * eps()
        @test RAS_DMFT.moment(Σ, 1) ≈ 0 atol = 1.0e-9
        # broadened
        Σ_IFG_lorentz = self_energy_IFG_lorentzian(Σ_H, C, W, δ)
        Σ_IFG_gauss = self_energy_IFG_gaussian(Σ_H, C, W, δ)
        @test Σ_IFG_lorentz != Σ_IFG_gauss
        @test iszero(imag(first(Σ_IFG_gauss))) # exponential decay results in zero
        @test minimum(imag(Σ_IFG_gauss)) < minimum(imag(Σ_IFG_lorentz)) # Gauss is steeper
    end # correlator
end # self-energy
