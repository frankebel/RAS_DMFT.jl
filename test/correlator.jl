using RAS_DMFT
using Fermions
using LinearAlgebra
using Test

@testset "correlator" begin
    # parameters
    n_bath = 101
    U = 2.0
    μ = U / 2
    ϵ_imp = -μ
    n_v_bit = 1
    n_c_bit = 1
    p = 1
    n_kryl = 100
    var = eps()

    Δ0 = hybridization_function_bethe_simple(n_bath)
    # Operators for positive frequencies. Negative ones are calculated by adjoint.
    fs = FockSpace(Orbitals(2 + n_v_bit + n_c_bit), FermionicSpin(1 // 2))
    c = annihilators(fs)
    n = occupations(fs)
    H_int = U * n[1, 1 // 2] * n[1, -1 // 2]
    d_dag = c[1, -1 // 2]' # d_↓^†
    q_dag = H_int * d_dag - d_dag * H_int  # q_↓^† = [H_int, d^†]
    O = [d_dag, q_dag]

    H, _, ψ0 = init_system(Δ0, H_int, ϵ_imp, n_v_bit, n_c_bit, p, var)

    @testset "Lanczos" begin
        # G+
        G_plus = correlator_plus(H, ψ0, d_dag, n_kryl)
        @test typeof(G_plus) === PolesSum{Float64, Float64}
        @test length(G_plus) === n_kryl
        @test all(>=(0), locations(G_plus))
        @test all(>=(0), weights(G_plus))
        @test moment(G_plus, 0) ≈ 0.5 atol = 2.0e2 * eps()
        # G-
        G_minus = correlator_minus(H, ψ0, d_dag', n_kryl)
        @test typeof(G_minus) === PolesSum{Float64, Float64}
        @test length(G_minus) === n_kryl
        @test all(<=(0), locations(G_minus))
        @test all(>=(0), weights(G_minus))
        @test moment(G_minus, 0) ≈ 0.5 atol = 2.0e2 * eps()

        # analytic moments at half-filling
        G = G_plus + G_minus
        @test moment(G, 0) ≈ 1 atol = 1.0e2 * eps()
        @test moment(G, 1) ≈ 0 atol = 1.0e3 * eps()
    end # Lanczos

    @testset "block Lanczos" begin
        # C+
        C_plus = correlator_plus(H, ψ0, O, n_kryl)
        @test typeof(C_plus) === PolesSumBlock{Float64, Float64}
        @test length(C_plus) == length(O) * n_kryl
        @test all(>=(0), locations(C_plus))
        # C-
        C_minus = correlator_minus(H, ψ0, map(adjoint, O), n_kryl)
        @test typeof(C_minus) === PolesSumBlock{Float64, Float64}
        @test length(C_minus) == length(O) * n_kryl
        @test all(<=(0), locations(C_minus))

        # analytic moments at half-filling
        C = transpose(C_minus) + C_plus
        m0 = [1 U / 2; U / 2 U^2 / 2]
        @test isapprox(moment(C, 0), m0; atol = 1.0e3 * eps())
        m1 = [0 U^2 / 4; U^2 / 4 U^3 / 4]
        @test isapprox(moment(C, 1), m1; atol = 1.0e7 * eps())

        # Hartree term
        O_H = O[1]' * O[2] + O[2] * O[1]'
        Σ_H = dot(ψ0, O_H, ψ0)
        @test Σ_H ≈ U / 2 rtol = 1.0e3 * eps()
    end # block Lanczos

    @testset "_warn_wrong_sign" begin
        # scalar
        P = PolesSum([-0.5, 1.0, 2.0], [2.0, 1.0, 1.0])
        @test_logs (:warn, r"C\+ has negative spectral weight 2\.0 on 1 pole\(s\)") RAS_DMFT._warn_wrong_sign(P, :plus)
        @test_logs (:warn, r"C\- has positive spectral weight 2\.0 on 2 pole\(s\)") RAS_DMFT._warn_wrong_sign(P, :minus)
        @test_nowarn RAS_DMFT._warn_wrong_sign(PolesSum([1.0, 2.0], [2.0, 1.0]), :plus)

        # block
        Pb = PolesSumBlock([-0.5, 1.0], [[1.0 0; 0 1], [2.0 0; 0 2]])
        @test_logs (:warn, r"C\+ has negative spectral weight 2\.0 on 1 pole\(s\)") RAS_DMFT._warn_wrong_sign(Pb, :plus)
        @test_logs (:warn, r"C\- has positive spectral weight 4\.0 on 1 pole\(s\)") RAS_DMFT._warn_wrong_sign(Pb, :minus)
        @test_throws ArgumentError RAS_DMFT._warn_wrong_sign(P, :foo)
    end # _warn_wrong_sign
end # correlator
