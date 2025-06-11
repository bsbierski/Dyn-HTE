#from SpinMC.jl 

"""
    InteractionMatrix

Represents a 3×3 matrix that defines the coupling between two spins.

# Fields
- `m11`, `m12`, `m13`: First row elements
- `m21`, `m22`, `m23`: Second row elements
- `m31`, `m32`, `m33`: Third row elements

Used to represent the exchange interaction tensor Jᵅᵝ where α,β ∈ {x,y,z}.
"""
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
"""
    UnitCell{D}

Represents the unit cell of a crystal lattice in D dimensions.

# Fields
- `primitive::NTuple{D,NTuple{D,Float64}}`: Primitive lattice vectors
- `basis::Vector{NTuple{D,Float64}}`: Positions of basis sites within the unit cell
- `interactions::Vector{Tuple{Int,Int,NTuple{D,Int},Matrix{Float64}}}`: Interactions between sites
- `interactionsOnsite::Vector{Matrix{Float64}}`: Onsite interactions for each basis site
- `interactionsField::Vector{Vector{Float64}}`: External field at each basis site

# Description
The unit cell serves as the building block for constructing the full lattice.
It contains information about the lattice vectors, basis sites, and interactions.
"""
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
"""
    Lattice{D,N}

Represents a D-dimensional lattice with N interactions per site.

# Fields
- `size::NTuple{D, Int}`: Linear extent of lattice in unit cells
- `length::Int`: Total number of sites
- `unitcell::UnitCell{D}`: Unit cell definition for the lattice
- `sitePositions::Vector{NTuple{D,Float64}}`: Positions of all sites
- `spins::Matrix{Float64}`: 3×N_sites matrix containing spin configuration
- `interactionSites::Vector{NTuple{N,Int}}`: For each site, indices of interacting sites
- `interactionMatrices::Vector{NTuple{N,InteractionMatrix}}`: Interaction matrices for each site
- `interactionOnsite::Vector{InteractionMatrix}`: Onsite interaction matrices
- `interactionField::Vector{NTuple{3,Float64}}`: Local field at each site

# Description
The Lattice struct represents the complete lattice with all its sites and interactions.
It is constructed by replicating the unit cell according to the specified size.
"""
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
"""
    Graph

Represents a connected weighted graph without isolated vertices.

# Fields
- `g::SimpleWeightedGraph{Int64, Int64}`: The underlying weighted graph

# Description
The Graph struct wraps a SimpleWeightedGraph, adding the constraint that
the graph must be connected. Used for representing vacuum graphs in
high-temperature series expansions.
"""
struct Graph ### connected vacuum graph
    g::SimpleWeightedGraph{Int64, Int64}
    function Graph(g::SimpleWeightedGraph{Int64, Int64})
        @assert is_connected(g)
        new(g)
    end

    function Graph(g::SimpleWeightedGraph{Int64, Float64})
        @assert is_connected(g)
        new(SimpleWeightedGraph(map(Int64, g.weights)))
    end
end

"""
    GraphG

Represents a connected graph with two external vertices (legs).

# Fields
- `g::SimpleWeightedGraph{Int64, Int64}`: The underlying weighted graph
- `jjp::Vector{Int64}`: Indices of the two external vertices [j,j']

# Description
GraphG represents graphs used for correlator calculations, with two
specified external vertices. The graph must be connected and the
external vertices must exist within the graph.
"""
mutable struct GraphG ### V-connected graphs for correlator G with TWO external legs at sites j,jp
    g::SimpleWeightedGraph{Int64, Int64}
    jjp::Vector{Int64}
    function GraphG(g::SimpleWeightedGraph{Int64, Int64},jjp)

        @assert is_connected(g)
        @assert nv(g) >= maximum(jjp)
        @assert length(jjp)==2

        new(g,jjp)
    end

    function GraphG(g::SimpleWeightedGraph{Int64, Float64},jjp)
        g = SimpleWeightedGraph(map(Int64, g.weights))
        @assert is_connected(g)
        @assert nv(g) >= maximum(jjp)
        @assert length(jjp)==2

        new(g,jjp)
    end
    
end


"""
    gG_properties

Stores properties of a GraphG object for efficient processing.

# Fields
- `order::Int`: Order of the graph (number of edges)
- `index::Int`: Index of the graph in its original collection
- `symmetry_factor::Int`: Symmetry factor of the graph
- `is_symmetric::Bool`: Whether the graph is symmetric under exchange of external vertices

# Description
This struct stores essential properties of a graph with external legs,
facilitating efficient computation without needing the full graph structure.
"""
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

"""
    unique_Graph

Represents a unique graph and its equivalent variants.

# Fields
- `ref_graph::GraphG`: Reference graph (representative)
- `distance::Int`: Graph distance between external vertices
- `gG_vec::Vector{gG_properties}`: Vector of properties of all multigraphs that are equivalent to the reference graph

# Description
Groups together a reference graph and information about other graphs
that are equivalent to it under symmetry operations.
"""
mutable struct unique_Graph
    ref_graph ::GraphG
    distance ::Int
    gG_vec ::Vector{gG_properties}

    function unique_Graph(data::Vector{Any})
        length(data) == 3 || throw(ArgumentError("Expected a vector of length 3"))
        return new(data[1], Int(data[3]),gG_properties.(data[2]))
    end

end

"""
    unique_Graphs

Collection of all unique_Graphs up to `max_order`.

# Fields
- `max_order::Int`: Maximum order of graphs in the collection
- `graphs::Vector{unique_Graph}`: Array of unique graphs with their equivalents

# Description
Organizes graphs by their unique topological structure, allowing
for efficient calculation by avoiding redundant computations.
"""
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

"""
    Dyn_HTE_Graphs

Main data structure for high-temperature series expansion calculations.

# Fields
- `S::Rational{Int}`: Spin length (e.g., 1/2 or 1)
- `unique_graphs::unique_Graphs`: Collection of unique graphs
- `c_dict::Vector{Vector{Vector{Rational{Int128}}}}`: Coefficients for each graph

# Description
Contains all the graph data and coefficients needed for high-temperature
series expansion calculations for a particular spin length.
"""
mutable struct Dyn_HTE_Graphs
    S :: Rational{Int}    #Spin-Length
    unique_graphs ::unique_Graphs #a dictionary ordering all graphs into equivalence classes 
    c_dict::Vector{Vector{Vector{Rational{Int128}}}} #a dictionary of all values "c" for all graphs 
end

"""
    Dyn_HTE_Lattice

Represents a lattice prepared for high-temperature series expansion calculations.

# Fields
- `name::String`: Name of the lattice geometry
- `lattice::Lattice`: The lattice structure
- `graph::SimpleGraph`: Graph representation of the lattice
- `basis_positions::Vector{<:Int}`: Indices of center/reference sites

# Description
Wraps a lattice with additional information needed for high-temperature
series expansion calculations, including the reference sites for correlations.
"""
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


