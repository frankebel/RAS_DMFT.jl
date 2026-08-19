using Fermions
using Fermions.Wavefunctions
using LinearAlgebra
using RAS_DMFT
using SparseArrays
using Test

@testset "natural orbitals" begin
    @testset "matrix transformation" begin
        Δ = hybridization_function_bethe_simple(11)
        m = arrowhead_matrix(Δ)
        H, n_occ = to_natural_orbitals(m)

        @test n_occ === 6
        @test size(H) === (12, 12)
        d = diag(H)
        @test norm(d[1:6] + d[7:12]) < 50 * eps() # PHS
        @test ishermitian(H)
        @test H ≈ [
            1.161007701055029e-16 0.18330573804110445 0.0 0.0 0.0 0.0 -0.42754884259276227 -0.18330573804110853 0.0 0.0 0.0 0.0
            0.18330573804110445 -0.5076151092043915 0.2355852387704299 0.0 0.0 0.0 0.1833057380411033 0.0 0.0 0.0 0.0 0.0
            0.0 0.2355852387704299 -0.5515345193381724 -0.2104776852190996 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 -0.2104776852190996 -0.6121205400242606 -0.17304126343374435 0.0 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 -0.17304126343374435 -0.7038886821335019 0.11528512946095112 0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.11528512946095112 -0.8454072119862861 0.0 0.0 0.0 0.0 0.0 0.0
            -0.42754884259276227 0.1833057380411033 0.0 0.0 0.0 0.0 5.457797046034841e-15 0.18330573804110967 0.0 0.0 0.0 0.0
            -0.18330573804110853 0.0 0.0 0.0 0.0 0.0 0.18330573804110967 0.5076151092043929 -0.2355852387704316 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.2355852387704316 0.551534519338173 0.21047768521910148 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.21047768521910148 0.612120540024263 -0.17304126343374962 0.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.17304126343374962 0.7038886821335045 -0.1152851294609553
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.1152851294609553 0.8454072119862802
        ] atol = 2.0e-13

        Δ = hybridization_function_bethe_simple(301)
        m = arrowhead_matrix(Δ)
        H, n_occ = to_natural_orbitals(m)
        @test n_occ === 151
        E = eigvals(H)
        @test norm(E[1:151] + reverse(E[152:end])) < sqrt(eps())

        # non-hermitian matrix
        m = rand(10, 10)
        @test_throws ArgumentError to_natural_orbitals(m)

        # complex matrix
        @test !hasmethod(to_natural_orbitals, Tuple{Matrix{ComplexF64}})

        # sites with zero hybridization, Löwdin must not return negative eigenvalues
        Δ = hybridization_function_bethe_grid(range(-2, 2; length = 31))
        to_natural_orbitals(arrowhead_matrix(Δ))
    end # matrix transformation

    @testset "operator" begin
        # H_nat1 has 0 in irrelevant entries, H_nat2 has all values set.
        # Both should still result in the same Hamiltonian.
        H_nat1 = [
            0 2 0 4 5 0
            2 8 9 10 0 0
            0 9 15 0 0 0
            4 10 0 22 23 0
            5 0 0 23 29 30
            0 0 0 0 30 36
        ]
        H_nat2 = [
            1 2 3 4 5 6
            2 8 9 10 11 12
            3 9 15 0 17 18
            4 10 0 22 23 24
            5 11 17 23 29 30
            6 12 18 24 30 36
        ]
        n_occ = 3
        U = 37
        μ = 38

        @testset "Operator" begin
            fs = FockSpace(Orbitals(12), FermionicSpin(1 // 2))
            n = occupations(fs)
            H_int = U * n[1, -1 // 2] * n[1, 1 // 2]
            H1 = natural_orbital_operator(H_nat1, H_int, -μ, fs, n_occ, 1, 1)
            H2 = natural_orbital_operator(H_nat2, H_int, -μ, fs, n_occ, 1, 1)
            @test H1 == H2

            c = annihilators(fs)
            n = occupations(fs)
            H3 =
                # impurity
                U * n[1, 1 // 2] * n[1, -1 // 2] +
                -38 * n[1, -1 // 2] +
                -38 * n[1, 1 // 2] +
                # mirror site
                22 * n[2, -1 // 2] +
                22 * n[2, 1 // 2] +
                # bit component
                8 * n[3, -1 // 2] +
                8 * n[3, 1 // 2] +
                29 * n[4, -1 // 2] +
                29 * n[4, 1 // 2] +
                # vector component
                15 * n[5, -1 // 2] +
                15 * n[5, 1 // 2] +
                36 * n[6, -1 // 2] +
                36 * n[6, 1 // 2] +
                # hopping i <-> b
                4 * c[1, -1 // 2]' * c[2, -1 // 2] +
                4 * c[2, -1 // 2]' * c[1, -1 // 2] +
                4 * c[1, 1 // 2]' * c[2, 1 // 2] +
                4 * c[2, 1 // 2]' * c[1, 1 // 2] +
                # hopping i <-> other
                2 * c[1, -1 // 2]' * c[3, -1 // 2] +
                2 * c[3, -1 // 2]' * c[1, -1 // 2] +
                2 * c[1, 1 // 2]' * c[3, 1 // 2] +
                2 * c[3, 1 // 2]' * c[1, 1 // 2] +
                5 * c[1, -1 // 2]' * c[4, -1 // 2] +
                5 * c[4, -1 // 2]' * c[1, -1 // 2] +
                5 * c[1, 1 // 2]' * c[4, 1 // 2] +
                5 * c[4, 1 // 2]' * c[1, 1 // 2] +
                # hopping b <-> other
                10 * c[2, -1 // 2]' * c[3, -1 // 2] +
                10 * c[3, -1 // 2]' * c[2, -1 // 2] +
                10 * c[2, 1 // 2]' * c[3, 1 // 2] +
                10 * c[3, 1 // 2]' * c[2, 1 // 2] +
                23 * c[2, -1 // 2]' * c[4, -1 // 2] +
                23 * c[4, -1 // 2]' * c[2, -1 // 2] +
                23 * c[2, 1 // 2]' * c[4, 1 // 2] +
                23 * c[4, 1 // 2]' * c[2, 1 // 2] +
                # hopping v
                9 * c[3, -1 // 2]' * c[5, -1 // 2] +
                9 * c[5, -1 // 2]' * c[3, -1 // 2] +
                9 * c[3, 1 // 2]' * c[5, 1 // 2] +
                9 * c[5, 1 // 2]' * c[3, 1 // 2] +
                # hopping c
                30 * c[4, -1 // 2]' * c[6, -1 // 2] +
                30 * c[6, -1 // 2]' * c[4, -1 // 2] +
                30 * c[4, 1 // 2]' * c[6, 1 // 2] +
                30 * c[6, 1 // 2]' * c[4, 1 // 2]
            @test H1 == H3

            # raise sites in bit component
            H1 = natural_orbital_operator(H_nat1, H_int, -μ, fs, n_occ, 2, 2)
            H2 = natural_orbital_operator(H_nat2, H_int, -μ, fs, n_occ, 2, 2)
            @test H1 == H2
            H3 =
                # impurity
                U * n[1, 1 // 2] * n[1, -1 // 2] +
                -38 * n[1, -1 // 2] +
                -38 * n[1, 1 // 2] +
                # mirror site
                22 * n[2, -1 // 2] +
                22 * n[2, 1 // 2] +
                # bit component
                8 * n[3, -1 // 2] +
                8 * n[3, 1 // 2] +
                15 * n[4, -1 // 2] +
                15 * n[4, 1 // 2] +
                29 * n[5, -1 // 2] +
                29 * n[5, 1 // 2] +
                36 * n[6, -1 // 2] +
                36 * n[6, 1 // 2] +
                # hopping i <-> b
                4 * c[1, -1 // 2]' * c[2, -1 // 2] +
                4 * c[2, -1 // 2]' * c[1, -1 // 2] +
                4 * c[1, 1 // 2]' * c[2, 1 // 2] +
                4 * c[2, 1 // 2]' * c[1, 1 // 2] +
                # hopping i <-> other
                2 * c[1, -1 // 2]' * c[3, -1 // 2] +
                2 * c[3, -1 // 2]' * c[1, -1 // 2] +
                2 * c[1, 1 // 2]' * c[3, 1 // 2] +
                2 * c[3, 1 // 2]' * c[1, 1 // 2] +
                5 * c[1, -1 // 2]' * c[5, -1 // 2] +
                5 * c[5, -1 // 2]' * c[1, -1 // 2] +
                5 * c[1, 1 // 2]' * c[5, 1 // 2] +
                5 * c[5, 1 // 2]' * c[1, 1 // 2] +
                # hopping b <-> other
                10 * c[2, -1 // 2]' * c[3, -1 // 2] +
                10 * c[3, -1 // 2]' * c[2, -1 // 2] +
                10 * c[2, 1 // 2]' * c[3, 1 // 2] +
                10 * c[3, 1 // 2]' * c[2, 1 // 2] +
                23 * c[2, -1 // 2]' * c[5, -1 // 2] +
                23 * c[5, -1 // 2]' * c[2, -1 // 2] +
                23 * c[2, 1 // 2]' * c[5, 1 // 2] +
                23 * c[5, 1 // 2]' * c[2, 1 // 2] +
                # # hopping v
                9 * c[3, -1 // 2]' * c[4, -1 // 2] +
                9 * c[4, -1 // 2]' * c[3, -1 // 2] +
                9 * c[3, 1 // 2]' * c[4, 1 // 2] +
                9 * c[4, 1 // 2]' * c[3, 1 // 2] +
                # # hopping c
                30 * c[5, -1 // 2]' * c[6, -1 // 2] +
                30 * c[6, -1 // 2]' * c[5, -1 // 2] +
                30 * c[5, 1 // 2]' * c[6, 1 // 2] +
                30 * c[6, 1 // 2]' * c[5, 1 // 2]
            @test H1 == H3

            # non-hermitian
            @test_throws ArgumentError natural_orbital_operator(
                rand(Int, 6, 6), H_int, -μ, fs, n_occ, 1, 1
            )
        end # Operator

        @testset "RASOperator" begin
            # some sites of each chain in bit component
            fs = FockSpace(Orbitals(4), FermionicSpin(1 // 2))
            n = occupations(fs)
            H_int = U * n[1, -1 // 2] * n[1, 1 // 2]
            H1 = natural_orbital_ras_operator(H_nat1, H_int, -μ, fs, n_occ, 1, 1, 0)
            H2 = natural_orbital_ras_operator(H_nat2, H_int, -μ, fs, n_occ, 1, 1, 0)
            @test typeof(H1) == typeof(H2)
            @test H1 == H2

            c = annihilators(fs)
            n = occupations(fs)
            H_bit =
                # impurity
                U * n[1, 1 // 2] * n[1, -1 // 2] +
                -38 * n[1, -1 // 2] +
                -38 * n[1, 1 // 2] +
                # mirror site
                22 * n[2, -1 // 2] +
                22 * n[2, 1 // 2] +
                # bit component
                8 * n[3, -1 // 2] +
                8 * n[3, 1 // 2] +
                29 * n[4, -1 // 2] +
                29 * n[4, 1 // 2] +
                # hopping i <-> b
                4 * c[1, -1 // 2]' * c[2, -1 // 2] +
                4 * c[2, -1 // 2]' * c[1, -1 // 2] +
                4 * c[1, 1 // 2]' * c[2, 1 // 2] +
                4 * c[2, 1 // 2]' * c[1, 1 // 2] +
                # hopping i <-> other
                2 * c[1, -1 // 2]' * c[3, -1 // 2] +
                2 * c[3, -1 // 2]' * c[1, -1 // 2] +
                2 * c[1, 1 // 2]' * c[3, 1 // 2] +
                2 * c[3, 1 // 2]' * c[1, 1 // 2] +
                5 * c[1, -1 // 2]' * c[4, -1 // 2] +
                5 * c[4, -1 // 2]' * c[1, -1 // 2] +
                5 * c[1, 1 // 2]' * c[4, 1 // 2] +
                5 * c[4, 1 // 2]' * c[1, 1 // 2] +
                # hopping b <-> other
                10 * c[2, -1 // 2]' * c[3, -1 // 2] +
                10 * c[3, -1 // 2]' * c[2, -1 // 2] +
                10 * c[2, 1 // 2]' * c[3, 1 // 2] +
                10 * c[3, 1 // 2]' * c[2, 1 // 2] +
                23 * c[2, -1 // 2]' * c[4, -1 // 2] +
                23 * c[4, -1 // 2]' * c[2, -1 // 2] +
                23 * c[2, 1 // 2]' * c[4, 1 // 2] +
                23 * c[4, 1 // 2]' * c[2, 1 // 2]
            H_mix = RASOperatorMixed{UInt64, Int}[]
            @test H1.opbit == H_bit
            @test H1.opmix == H_mix
            @test H1.zero == 2 * 15
            @test H1.one == SymTridiagonal(Int[], Int[])
            @test H1.two == sparse(Int[], Int[], Int[])
            @test H1.nbit == 4
            @test H1.nfilled == 1
            @test H1.nempty == 1
            @test H1.excitation == 0

            # single excitation
            a = repeat([30 - 15, 30 + 36], 2)
            b = zeros(Int, 3)
            one = SymTridiagonal(a, b)
            H1 = natural_orbital_ras_operator(H_nat1, H_int, -μ, fs, n_occ, 1, 1, 1)
            H2 = natural_orbital_ras_operator(H_nat2, H_int, -μ, fs, n_occ, 1, 1, 1)
            H_mix = [
                RASOperatorMixed(only((9 * c[3, -1 // 2]').terms), 0x0000000000000008, 1, 2),
                RASOperatorMixed(only((9 * c[3, -1 // 2]).terms), 0x0000000000000008, 2, 1),
                RASOperatorMixed(only((9 * c[3, 1 // 2]').terms), 0x0000000000000080, 1, 4),
                RASOperatorMixed(only((9 * c[3, 1 // 2]).terms), 0x0000000000000080, 4, 1),
                RASOperatorMixed(only((-30 * c[4, -1 // 2]').terms), 0x0000000000000000, 3, 1),
                RASOperatorMixed(only((-30 * c[4, -1 // 2]).terms), 0x0000000000000000, 1, 3),
                RASOperatorMixed(only((-30 * c[4, 1 // 2]').terms), 0x0000000000000000, 5, 1),
                RASOperatorMixed(only((-30 * c[4, 1 // 2]).terms), 0x0000000000000000, 1, 5),
            ]
            @test typeof(H1) == typeof(H2)
            @test H1 == H2
            @test H1.opbit == H_bit
            @test H1.opmix == H_mix
            @test H1.zero == 2 * 15
            @test H1.one == one
            @test H1.two == sparse(Int[], Int[], Int[])
            @test H1.nbit == 4
            @test H1.nfilled == 1
            @test H1.nempty == 1
            @test H1.excitation == 1

            # double excitation
            two = sparse(
                Diagonal(
                    [
                        30 - 15 + 36,
                        30 - 15 - 15,
                        30 - 15 + 36,
                        30 - 15 + 36,
                        30 + 36 + 36,
                        30 - 15 + 36,
                    ]
                ),
            )
            H1 = natural_orbital_ras_operator(H_nat1, H_int, -μ, fs, n_occ, 1, 1, 2)
            H2 = natural_orbital_ras_operator(H_nat2, H_int, -μ, fs, n_occ, 1, 1, 2)
            H_mix = [
                RASOperatorMixed(only((9 * c[3, -1 // 2]').terms), 0x0000000000000008, 1, 2),
                RASOperatorMixed(only((9 * c[3, -1 // 2]).terms), 0x0000000000000008, 2, 1),
                RASOperatorMixed(only((9 * c[3, 1 // 2]').terms), 0x0000000000000080, 1, 4),
                RASOperatorMixed(only((9 * c[3, 1 // 2]).terms), 0x0000000000000080, 4, 1),
                RASOperatorMixed(only((9 * c[3, 1 // 2]').terms), 0x0000000000000080, 2, 7),
                RASOperatorMixed(only((9 * c[3, 1 // 2]).terms), 0x0000000000000080, 7, 2),
                RASOperatorMixed(only((9 * c[3, -1 // 2]').terms), 0x0000000000000008, 3, 6),
                RASOperatorMixed(only((9 * c[3, -1 // 2]).terms), 0x0000000000000008, 6, 3),
                RASOperatorMixed(only((9 * c[3, 1 // 2]').terms), 0x0000000000000080, 3, 8),
                RASOperatorMixed(only((9 * c[3, 1 // 2]).terms), 0x0000000000000080, 8, 3),
                RASOperatorMixed(only((9 * c[3, -1 // 2]').terms), 0x0000000000000008, 4, 7),
                RASOperatorMixed(only((9 * c[3, -1 // 2]).terms), 0x0000000000000008, 7, 4),
                RASOperatorMixed(only((9 * c[3, -1 // 2]').terms), 0x0000000000000008, 5, 9),
                RASOperatorMixed(only((9 * c[3, -1 // 2]).terms), 0x0000000000000008, 9, 5),
                RASOperatorMixed(only((9 * c[3, 1 // 2]').terms), 0x0000000000000080, 5, 11),
                RASOperatorMixed(only((9 * c[3, 1 // 2]).terms), 0x0000000000000080, 11, 5),
                RASOperatorMixed(only((-30 * c[4, -1 // 2]').terms), 0x0000000000000000, 3, 1),
                RASOperatorMixed(only((-30 * c[4, -1 // 2]).terms), 0x0000000000000000, 1, 3),
                RASOperatorMixed(only((-30 * c[4, 1 // 2]').terms), 0x0000000000000000, 5, 1),
                RASOperatorMixed(only((-30 * c[4, 1 // 2]).terms), 0x0000000000000000, 1, 5),
                RASOperatorMixed(only((30 * c[4, -1 // 2]').terms), 0x0000000000000000, 6, 2),
                RASOperatorMixed(only((30 * c[4, -1 // 2]).terms), 0x0000000000000000, 2, 6),
                RASOperatorMixed(only((-30 * c[4, 1 // 2]').terms), 0x0000000000000000, 9, 2),
                RASOperatorMixed(only((-30 * c[4, 1 // 2]).terms), 0x0000000000000000, 2, 9),
                RASOperatorMixed(only((-30 * c[4, 1 // 2]').terms), 0x0000000000000000, 10, 3),
                RASOperatorMixed(only((-30 * c[4, 1 // 2]).terms), 0x0000000000000000, 3, 10),
                RASOperatorMixed(only((-30 * c[4, -1 // 2]').terms), 0x0000000000000000, 8, 4),
                RASOperatorMixed(only((-30 * c[4, -1 // 2]).terms), 0x0000000000000000, 4, 8),
                RASOperatorMixed(only((30 * c[4, 1 // 2]').terms), 0x0000000000000000, 11, 4),
                RASOperatorMixed(only((30 * c[4, 1 // 2]).terms), 0x0000000000000000, 4, 11),
                RASOperatorMixed(
                    only((-30 * c[4, -1 // 2]').terms), 0x0000000000000000, 10, 5
                ),
                RASOperatorMixed(only((-30 * c[4, -1 // 2]).terms), 0x0000000000000000, 5, 10),
            ]
            @test typeof(H1) == typeof(H2)
            @test H1 == H2
            @test H1.opbit == H_bit
            @test H1.opmix == H_mix
            @test H1.zero == 2 * 15
            @test H1.one == one
            @test H1.two == two
            @test H1.nbit == 4
            @test H1.nfilled == 1
            @test H1.nempty == 1
            @test H1.excitation == 2

            # more sites in bit component
            U1 = 4.0
            μ1 = U1 / 2
            fs = FockSpace(Orbitals(6), FermionicSpin(1 // 2))
            n = occupations(fs)
            H_int = U1 * n[1, -1 // 2] * n[1, 1 // 2]
            Δ = hybridization_function_bethe_simple(11)
            H_nat, n_occ1 = to_natural_orbitals(arrowhead_matrix(Δ))
            H = natural_orbital_ras_operator(H_nat, H_int, -μ1, fs, n_occ1, 2, 2, 2)
            @test length(H.opbit.terms) == 1 + 2 * 6 + 4 * 7
            @test length(H.opmix) == 96 # didn't calculate myself
            @test H.zero isa Float64
            @test size(H.one) == (2 * 6, 2 * 6)
            @test size(H.two) == (binomial(2 * 6, 2), binomial(2 * 6, 2))

            # no site of each chain in bit component
            # natural_orbital_ras_operator_zero
            fs = FockSpace(Orbitals(4), FermionicSpin(1 // 2)) # too many Orbitals on purpose
            n = occupations(fs)
            H_int = U * n[1, -1 // 2] * n[1, 1 // 2]
            H1 = RAS_DMFT._natural_orbital_ras_operator_zero(H_nat1, H_int, -μ, fs, n_occ, 2)
            H2 = RAS_DMFT._natural_orbital_ras_operator_zero(H_nat2, H_int, -μ, fs, n_occ, 2)
            @test typeof(H1) == typeof(H2)
            @test H1 == H2

            c = annihilators(fs)
            n = occupations(fs)
            H_bit =
                # impurity
                U * n[1, 1 // 2] * n[1, -1 // 2] +
                -38 * n[1, -1 // 2] +
                -38 * n[1, 1 // 2] +
                # mirror site
                22 * n[2, -1 // 2] +
                22 * n[2, 1 // 2] +
                # hopping i <-> b
                4 * c[1, -1 // 2]' * c[2, -1 // 2] +
                4 * c[2, -1 // 2]' * c[1, -1 // 2] +
                4 * c[1, 1 // 2]' * c[2, 1 // 2] +
                4 * c[2, 1 // 2]' * c[1, 1 // 2]
            zero = 2 * (8 + 15)
            a = repeat([-8, -15, 29, 36], 2) .+ zero
            b = [9, 0, 30, 0, 9, 0, 30]
            one = SymTridiagonal(a, b)

            @test H1.opbit == H_bit
            @test H1.zero === zero
            @test H1.one == one
            @test typeof(H1.two) === SparseMatrixCSC{Int, Int}
            @test size(H1.two) == (binomial(2 * 4, 2), binomial(2 * 4, 2))
            @test H1.nbit === 2
            @test H1.nfilled === 2
            @test H1.nempty === 2
            @test H1.excitation === 2
        end # RASOperator
    end # operator
end # natural orbitals
