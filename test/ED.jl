using RAS_DMFT
using RAS_DMFT.ED
using Fermions
using LinearAlgebra
using Test

@testset "ED" begin
    n_bath = 3
    U = 4.0
    ϵ_imp = -U / 2
    Δ = hybridization_function_bethe_simple(n_bath)
    fs = FockSpace(Orbitals(n_bath + 1), FermionicSpin(1 // 2))
    n = occupations(fs)
    H_int = U * n[1, 1 // 2] * n[1, -1 // 2]

    @testset "solve_impurity_ed" begin
        G = solve_impurity_ed(Δ, H_int, ϵ_imp)
        wgt = weights(G)
        n_sites = n_bath + 1
        @test length(G) ==
            2 * binomial(n_sites, n_sites ÷ 2) * binomial(n_sites, n_bath ÷ 2)
        @test moment(G, 0) ≈ 1 rtol = 10 * eps()
        # PHS
        @test norm(locations(G) + reverse(locations(G))) < 200 * eps()
        @test norm(wgt - reverse(wgt)) < 40 * eps()
    end # solve_impurity_ed
end # ED
