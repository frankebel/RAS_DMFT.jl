using RAS_DMFT
using Fermions
using LinearAlgebra
using Test

@testset "correlator" begin
    # parameters
    n_bath = 31
    U = 4.0
    μ = U / 2
    ϵ_imp = -μ
    n_v_bit = 1
    n_c_bit = 1
    e = 2
    n_kryl = 50
    var = eps()
    W = -10:0.002:10
    δ = 0.08

    Δ0 = hybridization_function_bethe_simple(n_bath)
    # Operators for positive frequencies. Negative ones are calculated by adjoint.
    fs = FockSpace(Orbitals(2 + n_v_bit + n_c_bit), FermionicSpin(1 // 2))
    c = annihilators(fs)
    n = occupations(fs)
    H_int = U * n[1, 1 // 2] * n[1, -1 // 2]
    d_dag = c[1, -1 // 2]' # d_↓^†
    q_dag = H_int * d_dag - d_dag * H_int  # q_↓^† = [H_int, d^†]
    O = [d_dag, q_dag]

    H, E0, ψ0 = init_system(Δ0, H_int, ϵ_imp, n_v_bit, n_c_bit, e, var)

    @testset "Lanczos" begin
        # G+
        G_plus = correlator_plus(H, ψ0, d_dag, n_kryl)
        @test typeof(G_plus) === PolesSum{Float64, Float64}
        @test length(G_plus) === 50
        @test all(>=(0), locations(G_plus))
        @test all(>=(0), weights(G_plus))
        @test moment(G_plus, 0) ≈ 0.5 atol = 100 * eps()
        # G-
        G_minus = correlator_minus(H, ψ0, d_dag', n_kryl)
        @test typeof(G_minus) === PolesSum{Float64, Float64}
        @test length(G_minus) === 50
        @test all(<=(0), locations(G_minus))
        @test all(>=(0), weights(G_minus))
        @test moment(G_minus, 0) ≈ 0.5 atol = 100 * eps()
        # symmetry: first moment must be zero
        G = G_plus + G_minus
        @test moment(G, 1) ≈ 0 atol = 100 * eps()
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

        # G±
        G_plus = PolesSum(copy(locations(C_plus)), map(i -> i[1, 1], weights(C_plus)))
        G_minus = PolesSum(copy(locations(C_minus)), map(i -> i[1, 1], weights(C_minus)))

        # half-filling
        @test moment(G_plus, 0) ≈ 0.5 rtol = 150 * eps()
        @test moment(G_minus, 0) ≈ 0.5 rtol = 150 * eps()

        # compare absolute moments of impurity Green's function
        m_pos = moments(G_plus, 0:10)
        m_neg = moments(G_minus, 0:10)
        ratio = m_pos ./ m_neg
        @test all(r -> isapprox(r, 1; atol = 500 * eps()), ratio[1:2:end])
        @test all(r -> isapprox(r, -1; atol = 500 * eps()), ratio[2:2:end])

        # Hartree term
        O_H = O[1]' * O[2] + O[2] * O[1]'
        Σ_H = dot(ψ0, O_H, ψ0)
        @test Σ_H ≈ U / 2 rtol = 100 * eps()

        # evaluation
        C = transpose(C_minus) + C_plus
        c = evaluate_lorentzian(C, W, δ)
        G = map(i -> i[1, 1], c)
        F = map(i -> i[1, 2], c)
        Σ = F ./ G
        @test all(i -> i <= 0, imag(Σ))
        @test real(Σ[cld(length(W), 2)]) ≈ U / 2 rtol = 500 * eps()
    end # block Lanczos
end # correlator
