module DynHTE

using Plots
using Graphs, SimpleWeightedGraphs, Parameters
using Symbolics, Polynomials, DifferentialEquations, LsqFit
using TaylorSeries, LinearAlgebra, SparseArrays, Combinatorics
using Random, JLD2
using PaddedViews, ToeplitzMatrices


import HypergeometricFunctions.pFq
import SpecialFunctions.gamma
import GraphRecipes.graphplot

# Add a function to access data paths
function data_dir()
    return joinpath(@__DIR__, "..", "data")
end

# Export the data directory function
export data_dir


# Include all component files
include("RobustPade.jl")
include("Structs.jl")
include("vf2_edited.jl")
include("Lattice.jl")
include("LatticeSymmetries.jl")
include("LatticeGraphs.jl")
include("GraphGeneration.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl")

# Export structs that users will need to interact with
export 
    # Main data structures
    Dyn_HTE_Graphs,
    Dyn_HTE_Lattice,
    Lattice,
    UnitCell,
    unique_Graphs,
    unique_Graph,
    GraphG,
    Graph,
    gG_properties

# Export the functions from ConvenienceFunctions.jl
export 
    # Spin string creation
    create_spin_string,
    
    # Matrix operations for correlators
    get_c_iipDyn_mat,
    get_c_iipEqualTime_mat,
    get_TGiip_Matsubara_xpoly,
    flipEvenIndexEntries,
    
    # Approximation methods
    get_pade,
    get_intDiffApprox,
    get_p_u,
    get_LinearTrafoToCoeffs_u,
    
    # Brillouin zone and k-space functions
    create_brillouin_zone_path,
    get_c_k,
    inverse_fourier_transform,
    get_c_k_subl,
    inverse_fourier_transform_subl,
    
    # Moment and continued fraction calculations
    get_moments_from_c_kDyn,
    fromMomentsToδ,
    contFrac,
    extrapolate_δvec,
    ContFracTerminator,
    get_extrapolation_params,
    JSwithTerminator,
    JS,
    get_JSkw_mat,
    extrapolate_series,
    find_divergence_point,
    
    # Additional required functions
    load_dyn_hte_graphs,
    getLattice,
    graphplot

end

