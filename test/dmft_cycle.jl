using RAS_DMFT
using Fermions
using LinearAlgebra
using Test

# Integration test: one complete step through DMFT.

@testset "DMFT cycle" begin
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
    W = range(-6, 6; length = 2001)
    tol_wgt = 1.0e-11

    Δ0 = hybridization_function_bethe_simple(n_bath)

    # Operators for positive frequencies. Negative ones are calculated by adjoint.
    fs = FockSpace(Orbitals(2 + n_v_bit + n_c_bit), FermionicSpin(1 // 2))
    c = annihilators(fs)
    n = occupations(fs)
    H_int = U * n[1, 1 // 2] * n[1, -1 // 2]
    d_dag = c[1, -1 // 2]' # d_↓^†
    q_dag = H_int * d_dag - d_dag * H_int  # q_↓^† = [H_int, d^†]

    # initialize system
    H, _, ψ0 = init_system(Δ0, H_int, ϵ_imp, n_v_bit, n_c_bit, e, var)

    # Hartree term and linear shift
    O_Σ_H = q_dag' * d_dag + d_dag * q_dag'
    Σ_H = dot(ψ0, O_Σ_H, ψ0)
    q_dag_tilde = q_dag - Σ_H * d_dag
    O = [q_dag_tilde, d_dag]

    # block correlator (impurity solver)
    C_plus = correlator_plus(H, ψ0, O, n_kryl)
    C_minus = correlator_minus(H, ψ0, map(adjoint, O), n_kryl)
    C = transpose(C_minus) + C_plus
    merge_small_weight!(C, tol_wgt)

    # dynamic part of self-energy
    Σ = PolesSum(self_energy_IFG(C), 1, 1)::PolesSum{Float64, Float64}
    merge_small_weight!(Σ, tol_wgt)
    @test moment(Σ, 0) ≈ U^2 / 4 atol = 2.0e2 * eps()
    @test moment(Σ, 1) ≈ 0 atol = 1.0e4 * eps()

    # new hybridization function on the grid `W`
    Δ = update_hybridization_function(Δ0, μ, Σ_H, Σ)
    merge_degenerate_poles!(Δ)
    merge_small_weight!(Δ, tol_wgt)
    Δ_new = to_grid(Δ, W)

    @test locations(Δ_new) == W
    @test moment(Δ_new, 0) ≈ 0.25 atol = 10 * eps()
    @test moment(Δ_new, 1) ≈ 0.0 atol = sqrt(eps())
end # DMFT cycle
