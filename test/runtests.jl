using DynHTE
using Test
using Graphs,SimpleWeightedGraphs,Polynomials,LinearAlgebra

# ContinuedFractionTerminator
import HypergeometricFunctions.pFq
import SpecialFunctions.gamma
import DynHTE: ContFracTerminator
#LatticeSymmetries
import DynHTE: sym_element, sym_group, translation_group, 
                is_element, mod!, neutral_elem, ∘, 
                add_element!, inverse, commutator, find_order,
                generate_closed_basis, generate_symmetry_group, 
                bond, flip_bond!, getSymmetryGroup, shiftRotation,floor_tol
#LatticeGraphs
import DynHTE: find_site_basis_label, find_graph_center,get_finite_Lattice, latticeToGraph

#Lattice.jl
import DynHTE: addBasisSite!, addInteraction!, setInteractionOnsite!, setField!, setSpin!, getSpin, getSitePosition

#Embedding 
import DynHTE: is_simple_isomorphic, e_fast, Calculate_Correlator_fast, SimpleGraph

@testset "DynHTE.jl" begin
    include("test_ConvenienceFunctions.jl")
    include("test_Embedding.jl")
    include("test_Structs.jl")
    include("test_Lattice.jl")
    include("test_LatticeGraphs.jl")
    include("test_LatticeSymmetries.jl")
end;
