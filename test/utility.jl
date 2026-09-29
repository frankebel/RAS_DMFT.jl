using Fermions
using Fermions.Lanczos
using Fermions.Wavefunctions
using LinearAlgebra
using RAS_DMFT
using Test

@testset "util" begin
    @testset "init_system" begin
        # parameters
        n_bath = 31
        U = 4.0
        μ = U / 2
        L_v = 1
        L_c = 1
        p = 2
        var = eps()

        E0_target = -21.527949990415227 # target ground state energy
        Δ = hybridization_function_bethe_simple(n_bath)
        fs = FockSpace(Orbitals(2 + L_v + L_c), FermionicSpin(1 // 2))
        n = occupations(fs)
        H_int = U * n[1, -1 // 2] * n[1, 1 // 2]
        H, E0, ψ0 = init_system(Δ, H_int, -μ, 0, L_v, L_c, p, var)
        Hψ = H * ψ0
        variance = Hψ ⋅ Hψ
        @test variance < var
        @test E0 ≈ E0_target rtol = 2.0e-15

        # no PHS by setting ϵ_mf ≠ 0
        H_p, E0_p, ψ_p = init_system(Δ, H_int, -μ, 0.5, L_v, L_c, p, var)
        Hψp = H_p * ψ_p
        variance_plus = Hψp ⋅ Hψp
        @test variance_plus < var
        H_m, E0_m, ψ_m = init_system(Δ, H_int, -μ, -0.5, L_v, L_c, p, var)
        Hψm = H_m * ψ_m
        variance_minus = Hψm ⋅ Hψm
        @test variance_minus < var
        @test E0_p ≈ E0_m rtol = 2.0e-15
        # shift destroys PHS
        @test !isapprox(E0_p, E0; rtol = 1.0e-9)

        # short Lanczos runs with many restarts reach the same ground state
        H_k, E0_k, ψ_k = init_system(
            Δ, H_int, -μ, 0, L_v, L_c, p, var; n_kryl = 5, n_max_restart = 100
        )
        Hψk = H_k * ψ_k
        @test Hψk ⋅ Hψk < var
        @test E0_k ≈ E0_target rtol = 2.0e-15
        # a single short restart cannot reach the variance
        @test_logs (:warn,) init_system(
            Δ, H_int, -μ, 0, L_v, L_c, p, var; n_kryl = 2, n_max_restart = 1
        )
    end # init system

    @testset "Kondo temperature" begin
        @test temperature_kondo(0.3, -0.1, 0.1) == 0.04297872341114842
        @test temperature_kondo(0.2, -0.1, 0.015) == 0.00020610334475146955
    end # Kondo temperature

    @testset "find_chemical_potential" begin
        n_tol = 1.0e-12

        # Each band is diagonal, so the bands decouple.
        function filling_exact(μ, H_ks, Σ_stat, locs, wgts; tol = 100 * eps())
            n_p = length(locs)
            total = 0.0
            for H_k in H_ks, i in eachindex(Σ_stat)
                A = zeros(1 + n_p, 1 + n_p)
                A[1, 1] = H_k[i] + Σ_stat[i] - μ
                for j in 1:n_p
                    A[1 + j, 1 + j] = locs[j]
                    A[1, 1 + j] = A[1 + j, 1] = sqrt(wgts[j][i])
                end
                F = eigen(Symmetric(A))
                for q in eachindex(F.values)
                    λ, wgt_q = F.values[q], abs2(F.vectors[1, q])
                    λ < -tol && (total += wgt_q)
                    abs(λ) <= tol && (total += 0.5 * wgt_q)
                end
            end
            return total / length(H_ks)
        end

        @testset "PHS" begin
            H_ks_diag = [[0.0]]
            Σ_stat_diag = [0.0]
            locs = [-1.0, 1.0]
            wgts_diag = [[0.5], [0.5]]
            H_ks = Diagonal.(H_ks_diag)
            Σ_stat = Diagonal(Σ_stat_diag)
            Σ_dyn = PolesSumBlock(locs, [Diagonal(wgt) for wgt in wgts_diag])

            μ1_ref = -0.5
            n1_ref = filling_exact(μ1_ref, H_ks_diag, Σ_stat_diag, locs, wgts_diag)
            μ1, n1 = find_chemical_potential(
                H_ks, Σ_stat, Σ_dyn, [1], n1_ref;
                μ_min = -1.0, μ_max = 1.0, μ_tol = 1.0e-13, b_max = 200, n_tol,
            )
            @test μ1 ≈ μ1_ref atol = 1.0e-12
            @test n1 ≈ n1_ref atol = 1.0e-12

            μ2_ref = 0.5
            n2_ref = filling_exact(μ2_ref, H_ks_diag, Σ_stat_diag, locs, wgts_diag)
            μ2, n2 = find_chemical_potential(
                H_ks, Σ_stat, Σ_dyn, [1], n2_ref;
                μ_min = -1.0, μ_max = 1.0, μ_tol = 1.0e-13, b_max = 200, n_tol,
            )
            @test μ2 ≈ μ2_ref atol = 1.0e-12
            @test n2 ≈ n2_ref atol = 1.0e-12

            @test n1 + n2 ≈ 1 atol = 1.0e-12 # PHS
            # filling has a jump at exactly zero
            # approach symmetric
            μ0, _ = find_chemical_potential(
                H_ks, Σ_stat, Σ_dyn, [1], 0.5;
                μ_min = -1.0, μ_max = 1.0, μ_tol = 1.0e-13, b_max = 200, n_tol,
            )
            @test μ0 ≈ 0 atol = 1.0e-12
            # approach from left
            μ0, n0 = find_chemical_potential(
                H_ks, Σ_stat, Σ_dyn, [1], 0.5;
                μ_min = -1.0, μ_max = 0.0, μ_tol = 1.0e-13, b_max = 200, n_tol,
            )
            @test μ0 ≈ 0 atol = 1.0e-12
            @test n0 ≈ 0.25 atol = 1.0e-12
            # approach from right
            μ0, n0 = find_chemical_potential(
                H_ks, Σ_stat, Σ_dyn, [1], 0.5;
                μ_min = 0.0, μ_max = 1.0, μ_tol = 1.0e-13, b_max = 200, n_tol,
            )
            @test μ0 ≈ 0 atol = 1.0e-12
            @test n0 ≈ 0.75 atol = 1.0e-12
        end # PHS

        @testset "basis and ordering" begin
            H_ks_diag = [[-0.7, 0.2, 1.1], [-0.4, 0.5, 0.9]]
            H_ks = Diagonal.(H_ks_diag)
            Σ_stat_diag = [0.15, -0.3, 0.4]
            Σ_stat = Diagonal(Σ_stat_diag)
            locs = [-1.2, 0.5]
            wgts = [[0.6, 0.0, 0.0], [0.0, 0.0, 0.35]]

            # dynamic self-energy only on bands 1 and 3, in either order
            idx13 = [1, 3]
            idx31 = [3, 1]
            Σ_dyn_13 = PolesSumBlock(locs, [Diagonal(wgt[idx13]) for wgt in wgts])
            Σ_dyn_31 = PolesSumBlock(locs, [Diagonal(wgt[idx31]) for wgt in wgts])

            # rotate by unitary matrix
            M = ComplexF64[1 2im -1; 0.5im -1 3; -2 1 0.5im]
            @test rank(M) == 3
            U = Matrix(qr(M).Q)
            H_rot = [Hermitian(U * Diagonal(H_k) * U') for H_k in H_ks_diag]
            Σ_stat_rot = Hermitian(U * Diagonal(Σ_stat_diag) * U')
            Σ_dyn_rot = PolesSumBlock(
                locs, [Hermitian(U * diagm(wgt) * U') for wgt in wgts]
            )

            # `Σ_dyn` has to match `idx`, and `idx` has to stay inside the bands.
            @test_throws DimensionMismatch find_chemical_potential(
                H_ks, Σ_stat, Σ_dyn_13, [1], 0.1
            )
            @test_throws ArgumentError find_chemical_potential(
                H_ks, Σ_stat, Σ_dyn_13, [1, 1], 0.1
            )
            @test_throws ArgumentError find_chemical_potential(
                H_ks, Σ_stat, Σ_dyn_13, [1, 4], 0.1
            )

            # The hard-coded values come from the diagonal case.
            # The same system in three layouts: the correlated block in either
            # order, and everything rotated into a dense complex basis.
            # All three have to return the chemical potential it came from.
            for (μ_target, n_fill) in (
                    (-0.5, 0.18015156588363956),
                    (0.75, 2.2394449917269812),
                    (1.7, 2.7786451479534993),
                )
                @test filling_exact(μ_target, H_ks_diag, Σ_stat_diag, locs, wgts) ≈
                    n_fill atol = 8 * eps()
                for (H, Σ_s, Σ_d, idx) in (
                        (H_ks, Σ_stat, Σ_dyn_13, idx13),
                        (H_ks, Σ_stat, Σ_dyn_31, idx31),
                        (H_rot, Σ_stat_rot, Σ_dyn_rot, 1:3),
                    )
                    μ, n = find_chemical_potential(
                        H, Σ_s, Σ_d, idx, n_fill;
                        μ_min = -4.0, μ_max = 6.0, μ_tol = 1.0e-13, b_max = 200, n_tol,
                    )
                    @test μ ≈ μ_target atol = 1.0e-12
                    @test n ≈ n_fill atol = 1.0e-12
                end
            end
        end # basis and ordering
    end # find_chemical_potential

    @testset "_issorted_and_unique" begin
        @test RAS_DMFT._issorted_and_unique(1:10)
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([1, 0]) # not sorted
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([0, 0, 1]) # not unique
        # duplicate zeros
        @test_throws ArgumentError RAS_DMFT._issorted_and_unique([-0.0, 0.0])
    end # _issorted_and_unique


end # util
