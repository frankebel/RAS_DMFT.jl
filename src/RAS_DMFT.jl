module RAS_DMFT

using Distributions
using Fermions
using Fermions.Bits
using Fermions.Lanczos
using Fermions.Wavefunctions
using HDF5
using LinearAlgebra
using SpecialFunctions
using Statistics

export
    # Types
    AbstractPoles,
    AbstractPolesContinuedFraction,
    AbstractPolesSum,
    PolesContinuedFraction,
    PolesContinuedFractionBlock,
    PolesSum,
    PolesSumBlock,

    # Functions
    CIWavefunction_singlet,
    Wavefunction_singlet,
    add_pole_at_zero!,
    amplitude,
    amplitudes,
    anderson_matrix,
    block_lanczos,
    block_lanczos_full_ortho,
    correlator,
    correlator_minus,
    correlator_plus,
    discretize_similar_weight,
    discretize_to_grid,
    equal_weight_discretization,
    evaluate_gaussian,
    evaluate_lorentzian,
    find_chemical_potential,
    flip_spectrum!,
    flip_spectrum,
    get_CI_parameters,
    greens_function_bethe_analytic,
    greens_function_bethe_equal_weight,
    greens_function_bethe_grid,
    greens_function_bethe_grid_hubbard3,
    greens_function_bethe_simple,
    greens_function_local,
    greens_function_partial,
    grid_interpolate,
    grid_log,
    ground_state!,
    hybridization_function_bethe_analytic,
    hybridization_function_bethe_equal_weight,
    hybridization_function_bethe_grid,
    hybridization_function_bethe_grid_hubbard3,
    hybridization_function_bethe_simple,
    imagKK,
    init_system,
    locations,
    mask_fe,
    merge_degenerate_poles!,
    merge_negative_locations_to_zero!,
    merge_negative_weight!,
    merge_small_weight!,
    moment,
    moments,
    natural_orbital_ci_operator,
    natural_orbital_operator,
    orgtr!,
    quasiparticle_weight,
    quasiparticle_weight_gaussian,
    read_hdf5,
    realKK,
    remove_zero_weight!,
    remove_zero_weight,
    scale,
    self_energy_IFG,
    self_energy_IFG_gaussian,
    self_energy_IFG_lorentzian,
    self_energy_dyson,
    slater_start,
    spectral_function_gauss,
    spectral_function_loggaussian,
    sytrd!,
    temperature_kondo,
    to_grid,
    to_natural_orbitals,
    update_hybridization_function,
    weight,
    weights,
    write_hdf5,
    δ_gaussian

include("bits.jl")
include("Poles/abstractpoles.jl")
include("Poles/abstractpolessum.jl")
include("Poles/polessum.jl")
include("Poles/polessumblock.jl")
include("Poles/abstractpolescontinuedfraction.jl")
include("Poles/polescontinuedfraction.jl")
include("Poles/polescontinuedfractionblock.jl")
include("Poles/conversion.jl")
include("io.jl")
include("sytrd.jl")
include("natural_orbitals.jl")
include("wavefunctions.jl")
include("kramers_kronig.jl")
include("grid.jl")
include("orthogonalization.jl")
include("utility.jl")
include("greens_function.jl")
include("hybridization_function.jl")
include("lanczos.jl")
include("block_lanczos.jl")
include("correlator.jl")
include("self_energy.jl")
include("update_hybridization_function.jl")
include("discretization.jl")
include("Combinatorics.jl")
include("Debug.jl")
include("ED.jl")

end
