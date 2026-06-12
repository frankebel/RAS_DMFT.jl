using RAS_DMFT
using Fermions
using Fermions.Wavefunctions
using LinearAlgebra
using Test

@testset "wavefunctions" begin
    @testset "Wavefunction_singlet" begin
        ψ = Wavefunction_singlet(Dict{UInt64, Float64}, 1, 2, 3, 4)
        d = Dict(
            UInt64(0b0000_111_00_1_01_0000_111_00_1_10) => 1 / sqrt(2),
            UInt64(0b0000_111_00_1_10_0000_111_00_1_01) => 1 / sqrt(2),
        )
        ϕ = Wavefunction(d)
        @test ψ == ϕ
    end # Wavefunction_singlet

    @testset "CIWavefunction_singlet" begin
        # excitation = 0
        ψ = CIWavefunction_singlet(Dict{UInt64, Float64}, 1, 2, 3, 4, 0)
        # vectors must be the equal but no egal
        foo = collect(values(ψ))
        @test foo[1] == foo[2]
        @test foo[1] !== foo[2]
        v = [1 / sqrt(2)]
        d = Dict(UInt64(0b00_1_01_00_1_10) => copy(v), UInt64(0b00_1_10_00_1_01) => copy(v))
        ϕ = CIWavefunction(d, 5, 3, 4, 0)
        @test ψ == ϕ

        # excitation = 1
        ψ = CIWavefunction_singlet(Dict{UInt64, Float64}, 1, 2, 3, 4, 1)
        v = zeros(1 + 2 * (3 + 4))
        v[1] = 1 / sqrt(2)
        d = Dict(UInt64(0b00_1_01_00_1_10) => copy(v), UInt64(0b00_1_10_00_1_01) => copy(v))
        ϕ = CIWavefunction(d, 5, 3, 4, 1)
        @test ψ == ϕ

        # excitation = 2
        ψ = CIWavefunction_singlet(Dict{UInt64, Float64}, 1, 2, 3, 4, 2)
        v = zeros(1 + 14 + 14 * 13 ÷ 2)
        v[1] = 1 / sqrt(2)
        d = Dict(UInt64(0b00_1_01_00_1_10) => copy(v), UInt64(0b00_1_10_00_1_01) => copy(v))
        ϕ = CIWavefunction(d, 5, 3, 4, 2)
        @test ψ == ϕ
    end # CIWavefunction_singlet

    @testset "ground state" begin
        # parameters
        n_bath = 31
        U = 4.0
        μ = U / 2
        L_v = 1
        L_c = 1
        p = 2
        n_kryl = 5
        n_sites = 1 + n_bath

        Δ = hybridization_function_bethe_simple(n_bath)
        H_nat, n_occ = to_natural_orbitals(Array(Δ))
        n_bit, V_v, V_c = get_CI_parameters(n_sites, n_occ, L_c, L_v)

        # CIWavefunction
        fs = FockSpace(Orbitals(n_bit), FermionicSpin(1 // 2))
        n = occupations(fs)
        H_int = U * n[1, -1 // 2] * n[1, 1 // 2]
        H = natural_orbital_ci_operator(H_nat, H_int, -μ, fs, n_occ, L_v, L_c, p)
        ψ_start = CIWavefunction_singlet(Dict{UInt64, Float64}, L_v, L_c, V_v, V_c, p)

        # 10 steps total
        E0, ψ0 = RAS_DMFT.ground_state!(H, ψ_start, n_kryl, 10, 0)
        @test E0 ≈ -21.52794995443462 rtol = 2.0e-13
        foo = H * ψ0
        var = foo ⋅ foo
        @test var < 3.0e-8

        # 20 steps total
        E0, ψ0 = RAS_DMFT.ground_state!(H, ψ0, n_kryl, 10, 0)
        @test E0 ≈ -21.52794999041698 rtol = 2.0e-13
        foo = H * ψ0
        var = foo ⋅ foo
        @test var < 3.0e-13

        # 30 steps total
        E0, ψ0 = RAS_DMFT.ground_state!(H, ψ0, n_kryl, 10, 0)
        @test E0 ≈ -21.52794999041725 rtol = 2.0e-13
        foo = H * ψ0
        var = foo ⋅ foo
        @test var < eps()

        # calculate full variance in Wavefunction without L,p approximation
        fs = FockSpace(Orbitals(n_sites), FermionicSpin(1 // 2))
        n = occupations(fs)
        H_int = U * n[1, -1 // 2] * n[1, 1 // 2]
        H_wf = natural_orbital_operator(H_nat, H_int, -μ, fs, n_occ, L_v, L_c)
        Fermions.shift_spectrum!(H_wf, E0)
        ψ0_wf = Wavefunction(ψ0)
        foo = H_wf * ψ0_wf
        var = foo ⋅ foo
        @test var < 2.0e-4
    end # ground state
end # wavefunctions
