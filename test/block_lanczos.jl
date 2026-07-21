using RAS_DMFT
using Fermions
using Fermions.Wavefunctions
using LinearAlgebra
using Test

@testset "block_lanczos" begin
    @testset "block_lanczos" begin
        # parameters
        n_bath = 11
        U = 4.0
        μ = U / 2
        n_v_bit = 1
        n_c_bit = 1
        e = 2
        n_kryl = 5
        var = eps()
        # initial system
        Δ = hybridization_function_bethe_simple(n_bath, 2)
        fs = FockSpace(Orbitals(2 + n_v_bit + n_c_bit), FermionicSpin(1 // 2))
        c = annihilators(fs)
        n = occupations(fs)
        H_int = U * n[1, 1 // 2] * n[1, -1 // 2]
        d_dag = c[1, -1 // 2]' # d_↓^†
        q_dag = H_int * d_dag - d_dag * H_int  # q_↓^† = [H_int, d^†]

        H, _, ψ0 = init_system(Δ, H_int, -μ, n_v_bit, n_c_bit, e, var)
        v1 = d_dag * ψ0
        v2 = q_dag * ψ0
        V0 = [v1 v2]
        # Löwdin orthonormalization
        W, _ = RAS_DMFT._orthonormalize_SVD(V0)
        # Block Lanczos
        a, b = RAS_DMFT.block_lanczos(H, W, n_kryl)
        @test length(a) == n_kryl
        @test length(b) == n_kryl - 1
        @test all(ishermitian, a)
        @test all(ishermitian, b)

        X = zeros(2 * n_kryl, 2 * n_kryl)
        for i in 1:(n_kryl - 1)
            X[(2 * i - 1):(2 * i), (2 * i - 1 + 2):(2 * i + 2)] = b[i]' # upper diagonal
            X[(2 * i - 1):(2 * i), (2 * i - 1):(2 * i)] = a[i] # main diagonal
            X[(2 * i + 1):(2 * i + 2), (2 * i - 1):(2 * i)] = b[i] # lower diagonal
        end
        X[(end - 1):end, (end - 1):end] = a[end] # last element
        E = eigvals(X)
        E_ref = [
            0.214059965436778
            0.5971189953214974
            0.9281541797197747
            1.132106853579721
            1.6228297233186255
            1.8127643102174924
            2.4744861501689104
            3.354246572428307
            4.674799177427353
            6.074676075106131
        ]
        @test norm(E - E_ref) < 3.0e3 * eps()
    end # block_lanczos

    @testset "block_lanczos_full_ortho" begin
        H = Diagonal([1, 1, 2, 2])
        Q1 = [1 0; 0 1 / sqrt(2); 0 0; 0 1 / sqrt(2)]
        A, B, Q = block_lanczos_full_ortho(H, Q1, 4)
        @test length(A) == 2
        @test norm(A[1] - [1 0; 0 1.5]) < 10 * eps()
        @test norm(A[2] - [0 0; 0 1.5]) < 10 * eps()
        @test length(B) == 1
        @test norm(B[1] - [0 0; 0 0.5]) < 10 * eps()
        @test length(Q) == 2
        @test Q[1] == Q1
        @test norm(Q[2] - [0 0; 0 -1 / sqrt(2); 0 0; 0 1 / sqrt(2)]) < 10 * eps()
    end # block_lanczos_full_ortho
end # block_lanczos
