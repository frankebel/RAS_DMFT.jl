module ED

using ..RAS_DMFT
using Fermions
using LinearAlgebra
using SparseArrays

export
    # Functions
    solve_impurity_ed

"""
    solve_impurity_ed(Δ::PolesSum, H_int::Operator, ϵ_imp::Real)

Solve AIM in exact diagonalization.
"""
function solve_impurity_ed(Δ::PolesSum, H_int::Operator, ϵ_imp::Real)
    # Get Hamiltonian
    H_nat, n_occ = to_natural_orbitals(arrowhead_matrix(Δ))
    n_sites = size(H_nat, 1)
    fs = FockSpace(Orbitals(n_sites), FermionicSpin(1 // 2))
    H = natural_orbital_operator(
        H_nat, H_int, ϵ_imp, fs, n_occ, n_occ - 1, n_sites - n_occ - 1
    )

    c = annihilators(fs)
    qn = NSzSet(fs)

    # half-filling
    block = Block(qn, (n_occ, n_occ)) # half-filling
    h = Array(BlockOper(H, block))
    @time "decomposition half-filling" E0, V0 = eigen(h)
    ψ0 = V0[:, 1] # ground state
    e0 = E0[1]

    # block with one less particle
    foo = BlockOper(c[1, -1 // 2], block) * ψ0
    h_minus = Array(BlockOper(H, Block(qn, (n_occ - 1, n_occ)))) # 1 less in spin ↓
    @time "decomposition H minus" a_minus, V_minus = eigen(h_minus)
    b_minus = [dot(view(V_minus, :, i), foo) for i in axes(V_minus, 2)]
    b_minus .= abs2.(b_minus)
    a_minus .-= e0
    a_minus .*= -1
    # sort increasing
    reverse!(a_minus)
    reverse!(b_minus)

    # block with one more particle
    foo = BlockOper(c[1, -1 // 2]', block) * ψ0
    h_plus = Array(BlockOper(H, Block(qn, (n_occ + 1, n_occ)))) # 1 more in spin ↓
    @time "decomposition H plus" a_plus, V_plus = eigen(h_plus)
    b_plus = [dot(view(V_plus, :, i), foo) for i in axes(V_plus, 2)]
    b_plus .= abs2.(b_plus)
    a_plus .-= e0

    # impurity GF
    return PolesSum([a_minus; a_plus], [b_minus; b_plus])
end

end
