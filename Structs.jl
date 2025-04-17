using Graphs, SimpleWeightedGraphs, Parameters

#from SpinMC.jl 
struct InteractionMatrix
    m11::Float64
    m12::Float64
    m13::Float64
    m21::Float64
    m22::Float64
    m23::Float64
    m31::Float64
    m32::Float64
    m33::Float64
end

#### UnitCell.jl
struct UnitCell{D}
    primitive::NTuple{D,NTuple{D,Float64}}
    basis::Vector{NTuple{D,Float64}}
    interactions::Vector{Tuple{Int,Int,NTuple{D,Int},Matrix{Float64}}} #interactions specified as (basis1,basis2,offsetPrimitive,M)
    interactionsOnsite::Vector{Matrix{Float64}}
    interactionsField::Vector{Vector{Float64}}

    UnitCell(a1::NTuple{1,Float64}) = new{1}((a1,), Vector{NTuple{1,Float64}}(undef,0), Vector{Tuple{Int,Int,NTuple{1,Int},Matrix{Float64}}}(undef,0), Vector{Matrix{Float64}}(undef,0), Vector{Vector{Float64}}(undef,0))
    UnitCell(a1::NTuple{2,Float64}, a2::NTuple{2,Float64}) = new{2}((a1,a2), Vector{NTuple{2,Float64}}(undef,0), Vector{Tuple{Int,Int,NTuple{2,Int},Matrix{Float64}}}(undef,0), Vector{Matrix{Float64}}(undef,0), Vector{Vector{Float64}}(undef,0))
    UnitCell(a1::NTuple{3,Float64}, a2::NTuple{3,Float64}, a3::NTuple{3,Float64}) = new{3}((a1,a2,a3), Vector{NTuple{3,Float64}}(undef,0), Vector{Tuple{Int,Int,NTuple{3,Int},Matrix{Float64}}}(undef,0), Vector{Matrix{Float64}}(undef,0), Vector{Vector{Float64}}(undef,0))
    UnitCell(primitives...) = new{length(primitives)}(primitives, Vector{NTuple{length(primitives),Float64}}(undef,0), Vector{Tuple{Int,Int,NTuple{length(primitives),Int},Matrix{Float64}}}(undef,0), Vector{Matrix{Float64}}(undef,0), Vector{Vector{Float64}}(undef,0))
end

#### Lattice.jl
mutable struct Lattice{D,N}
    size::NTuple{D, Int} #linear extent of the lattice in number of unit cells
    length::Int #Number of sites N_sites
    unitcell::UnitCell{D}
    sitePositions::Vector{NTuple{D,Float64}}

    spins::Matrix{Float64} #3*N_sites matrix containing the spin configuration

    interactionSites::Vector{NTuple{N,Int}} #list of length N_sites, for every site contains all interacting sites
    interactionMatrices::Vector{NTuple{N,InteractionMatrix}} #list of length N_sites, for every site contains all interaction matrices
    interactionOnsite::Vector{InteractionMatrix} #list of length N_sites, for every site contains the local onsite interaction matrix
    interactionField::Vector{NTuple{3,Float64}} #list of length N_sites, for every site contains the local field
    Lattice(D,N) = new{D,N}()
end


###### graph structures, all Graphs are: connected (and thus free of vertices without any edges)
struct Graph ### connected vacuum graph
    g::SimpleWeightedGraph{Int64, Int64}
    function Graph(g)
        @assert is_connected(g)
        new(g)
    end
end

mutable struct GraphG ### V-connected graphs for correlator G with TWO external legs at sites j,jp
    g::SimpleWeightedGraph{Int64, Int64}
    jjp::Vector{Int64}
    function GraphG(g,jjp)

        @assert is_connected(g)
        @assert nv(g) >= maximum(jjp)
        @assert length(jjp)==2

        new(g,jjp)
    end
end

mutable struct gG_properties
    order::Int
    index::Int
    symmetry_factor::Int
    is_symmetric::Bool

    function gG_properties(data::Vector{Int64})
        length(data) == 4 || throw(ArgumentError("Expected a vector of length 4"))
        return new(Int(data[1]), Int(data[2]), Int(data[3]),Bool(data[4]))
    end
end

mutable struct unique_Graph
    ref_graph ::GraphG
    distance ::Int
    gG_vec ::Vector{gG_properties}

    function unique_Graph(data::Vector{Any})
        length(data) == 3 || throw(ArgumentError("Expected a vector of length 3"))
        return new(data[1], Int(data[3]),gG_properties.(data[2]))
    end

end

mutable struct unique_Graphs
    max_order ::Int
    graphs ::Vector{unique_Graph}

    function unique_Graphs(data::Vector{Any})
        length(data) == 2 || throw(ArgumentError("Expected a vector of length 2"))
        return new(data[1], unique_Graph.(data[2]))
    end

    function unique_Graphs(max_order,graphs)
        return new(max_order, graphs)
    end

end

mutable struct Dyn_HTE_Graphs
    S :: Rational{Int}    #Spin-Length
    unique_graphs ::unique_Graphs #a dictionary ordering all graphs into equivalence classes 
    c_dict::Vector{Vector{Vector{Rational{Int128}}}} #a dictionary of all values "c" for all graphs 
end

mutable struct Dyn_HTE_Lattice
    name::String
    lattice::Lattice
    graph::SimpleGraph
    basis_positions::Vector{<:Int}
end

mutable struct unique_Graph_precalc
    ref_graph ::GraphG
    distance ::Int
    graph_value ::Matrix{Rational{Int128}}
    
    function unique_Graph_precalc(ref_graph,distance,graph_value)
        return new(ref_graph,distance,graph_value)
    end

end

mutable struct unique_Graphs_precalc
    max_order ::Int
    graphs ::Vector{unique_Graph_precalc}

    function unique_Graphs_precalc(max_order,graphs)
        return new(max_order, graphs)
    end

end


