import Base.-

"""
    latticeToGraph(lattice::Lattice) -> SimpleGraph{Int}

Convert a Lattice struct into a SimpleGraph representation.

# Arguments
- `lattice::Lattice`: The lattice structure to convert

# Returns
- `SimpleGraph{Int}`: Graph where vertices correspond to lattice sites and edges to interactions

# Description
Creates a graph representation of the lattice where each vertex corresponds to a site
in the lattice, and edges connect sites that interact with each other according to
the `interactionSites` field of the lattice.

# Examples
```julia
graph = latticeToGraph(lattice)
```
"""
function latticeToGraph(lattice::Lattice)::SimpleGraph{Int}
   
    g = SimpleGraph(lattice.length)
    ints = lattice.interactionSites
    for v in 1:nv(g)
        z = length(ints[v])
        for bond = 1:z
            if ints[v][bond] != v
                add_edge!(g, v, ints[v][bond])
            end
        end
    end
    
    return g
end


"""
    get_finite_Lattice(L::Int, geometry::String; PBC::Bool = true) -> (Lattice, SimpleGraph)

Create a lattice and corresponding graph with specified geometry and size.
Can be modified to include custom lattice geometries.

# Arguments
- `L::Int`: Linear size of the lattice
- `geometry::String`: Type of lattice geometry
- `PBC::Bool = true`: Whether to use periodic boundary conditions

# Returns
- `Tuple{Lattice, SimpleGraph}`: A tuple containing the lattice and its graph representation

# Supported Geometries
- "chain": 1D chain
- "square": 2D square lattice
- "simple_cubic": 3D cubic lattice
- "triang": 2D triangular lattice
- "honeycomb": 2D honeycomb lattice
- "kagome": 2D kagome lattice
- "pyrochlore": 3D pyrochlore lattice

# Examples
```julia
# Create a 10×10 square lattice with periodic boundaries
lattice, graph = get_finite_Lattice(10, "square")

# Create a honeycomb lattice without periodic boundaries
lattice, graph = get_finite_Lattice(8, "honeycomb", PBC=false)
```
"""
function get_finite_Lattice(L::Int,geometry::String; PBC::Bool = true)
  
    if geometry == "chain" ### chain lattice
        a1 = (1, 0)
        a2 = (0, 1)

        uc = UnitCell(a1,a2)
        b0 = addBasisSite!(uc, (0.0,0.0))

        addInteraction!(uc, b0, b0, (1,0))

        l = (L,1)

    elseif geometry == "square" ### Square lattice
        a1 = (1, 0)
        a2 = (0, 1)
        uc = UnitCell(a1,a2)

        b0 = addBasisSite!(uc, (0.0, 0.0))

        addInteraction!(uc, b0, b0, (1, 0))
        addInteraction!(uc, b0, b0, (0, 1))

        l = (L, L)

    elseif geometry == "simple_cubic" ### Square lattice
        a1 = (1, 0, 0)
        a2 = (0, 1 ,0)
        a3 = (0, 0, 1)
        uc = UnitCell(a1,a2,a3)

        b0 = addBasisSite!(uc, (0.0, 0.0, 0.0))

        addInteraction!(uc, b0, b0, (1, 0, 0))
        addInteraction!(uc, b0, b0, (0, 1 ,0))
        addInteraction!(uc, b0, b0, (0, 0, 1))

        l = (L, L ,L)

        
    elseif geometry == "triang"  ### Triangular lattice
        a1 = (1/2, sqrt(3)/2)
        a2 = (1/2, -sqrt(3)/2)
        uc = UnitCell(a1,a2)

        b0 = addBasisSite!(uc, (0.0, 0.0))

        addInteraction!(uc, b0, b0, (1, 0))
        addInteraction!(uc, b0, b0, (0, 1))
        addInteraction!(uc, b0, b0, (1, 1))

        l = (L, L)
 
    elseif geometry == "honeycomb" ### Honeycomb lattice
        a1 = (3/2, sqrt(3)/2)
        a2 = (3/2, -sqrt(3)/2)
        uc = UnitCell(a1,a2)

        b1 = addBasisSite!(uc, (0.0, 0.0))
        b2 = addBasisSite!(uc, (1.0, 0.0))

        addInteraction!(uc, b1, b2, (0, 0))
        addInteraction!(uc, b1, b2, (0, -1))
        addInteraction!(uc, b1, b2, (-1, 0))

        l = (L, L)
        
    elseif geometry == "kagome" ### Honeycomb lattice
        a1 = (1, sqrt(3))
        a2 = (1, -sqrt(3))
        uc = UnitCell(a1,a2)

        b1 = addBasisSite!(uc, (0.0, 0.0))
        b2 = addBasisSite!(uc, -1/2 .* a1)
        b3 = addBasisSite!(uc, 1/2 .* a2)

        addInteraction!(uc, b1, b2, (0, 0))
        addInteraction!(uc, b1, b3, (0, 0))
        addInteraction!(uc, b1, b2, (1, 0))
        addInteraction!(uc, b1, b3, (0, -1))

        addInteraction!(uc, b3, b2, (0, 0))
        addInteraction!(uc, b3, b2, (1, 1))

        

        l = (L, L)

    elseif geometry == "pyrochlore" ### Pyrochlore lattice
        a1 = (0,1/2,1/2)
        a2 = (1/2,0,1/2)
        a3 = (1/2,1/2,0)
        uc = UnitCell(a1,a2,a3)

        b0 = addBasisSite!(uc, (0.0, 0.0, 0.0))
        b1 = addBasisSite!(uc, (0.0, 1/4,1/4))
        b2 = addBasisSite!(uc, (1/4, 0.0 ,1/4))
        b3 = addBasisSite!(uc, (1/4,1/4, 0.0))

        addInteraction!(uc, b0, b1, (0, 0, 0))
        addInteraction!(uc, b0, b2, (0, 0, 0))
        addInteraction!(uc, b0, b3, (0, 0, 0))
        addInteraction!(uc, b0, b1, (-1, 0, 0))
        addInteraction!(uc, b0, b2, (0, -1, 0))
        addInteraction!(uc, b0, b3, (0, 0, -1))

        addInteraction!(uc, b1, b2, (0, 0, 0))
        addInteraction!(uc, b1, b3, (0, 0, 0))
        addInteraction!(uc, b2, b3, (0, 0, 0))
        addInteraction!(uc, b1, b2, (1, -1, 0))
        addInteraction!(uc, b1, b3, (1, 0, -1))
        addInteraction!(uc, b2, b3, (0, 1, -1))

        l = (L, L, L)

    else 
        error("geometry: " * geometry * " not yet implemented") 
    end

    ### create lattice with or without periodic Boundary conditions
    if PBC 
        lattice = Lattice(uc, l);
    else
        lattice = LatticeNoPBC(uc, l);
    end

    ### create graph for the lattice and return
    graph = latticeToGraph(lattice)
    return (lattice,graph)
end

"""
    -(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64}) -> Tuple{Float64, Float64}

Subtract two 2D coordinate tuples element-wise.

# Arguments
- `a::Tuple{Float64, Float64}`: First coordinate
- `b::Tuple{Float64, Float64}`: Second coordinate

# Returns
- `Tuple{Float64, Float64}`: Result of a - b

# Examples
```julia
diff = (3.0, 4.0) - (1.0, 2.0)  # Returns (2.0, 2.0)
```
"""
function -(a::Tuple{Float64, Float64},b::Tuple{Float64, Float64})::Tuple{Float64, Float64}
   
    return (a[1]-b[1],a[2]-b[2])
end

"""
    -(a::Tuple{Float64, Float64, Float64}, b::Tuple{Float64, Float64, Float64}) -> Tuple{Float64, Float64, Float64}

Subtract two 3D coordinate tuples element-wise.

# Arguments
- `a::Tuple{Float64, Float64, Float64}`: First coordinate
- `b::Tuple{Float64, Float64, Float64}`: Second coordinate

# Returns
- `Tuple{Float64, Float64, Float64}`: Result of a - b

# Examples
```julia
diff = (3.0, 4.0, 5.0) - (1.0, 2.0, 3.0)  # Returns (2.0, 2.0, 2.0)
```
"""
function -(a::Tuple{Float64, Float64, Float64},b::Tuple{Float64, Float64, Float64})::Tuple{Float64, Float64, Float64}
    
    return (a[1]-b[1],a[2]-b[2],a[3]-b[3])
end


"""
    find_graph_center(graph) -> Vector{Int}

Find vertices at the center of a graph based on average shortest path distances.

# Arguments
- `graph`: A graph structure supporting Graphs.floyd_warshall_shortest_paths

# Returns
- `Vector{Int}`: Indices of vertices with the minimum average distance to all other vertices

# Description
Computes the average shortest path distance from each vertex to all others,
then returns the vertices with the minimum average distance. These vertices
are considered the "center" of the graph in terms of connectivity.

# Examples
```julia
# Find the central vertices in a lattice graph
centers = find_graph_center(lattice_graph)
```
"""
function find_graph_center(graph)
    

    n = nv(graph)  # Number of vertices
    distances = Graphs.floyd_warshall_shortest_paths(graph).dists  # All-pairs shortest path distances
    
    function mean(x)
        sum(x)/length(x)
    end

    # Compute the average distance for each vertex
    avg_distances = [mean(filter(x -> x < typemax(Int), distances[i, :])) for i in 1:n]
    
    # Find the minimum average distance
    min_avg_distance = minimum(avg_distances)
    
    # Find the vertices with the minimum average distance
    min_vertices = findall(x -> x == min_avg_distance, avg_distances)
    
    return min_vertices
end


"""
    getLattice(L::Int, geometry::String) -> Dyn_HTE_Lattice

Generate a lattice where all sites are at most L steps away from center sites.

# Arguments
- `L::Int`: Maximum distance from center sites
- `geometry::String`: Type of lattice geometry (e.g., "chain", "square", "honeycomb")

# Returns
- `Dyn_HTE_Lattice`: A structure containing the modified lattice, its graph 
  representation, and identified center sites. See [`Dyn_HTE_Lattice`](@ref).

# Description
Creates a finite lattice centered around a reference site(s), including only
sites that are within L steps (graph distance) from the center. For a chain
lattice, a shortcut implementation is used with the middle site as center.

Uses [`get_finite_Lattice`](@ref) to generate the initial lattice.
To implement additional geometries, modify the `get_finite_Lattice` function.

# Examples
```julia
# Get a square lattice with sites at most 5 steps from center
lattice_ball = getLattice(5, "square")

# Get a honeycomb lattice with sites at most 4 steps from center
honeycomb_ball = getLattice(4, "honeycomb")
```
"""
function getLattice(L::Int,geometry::String)::Dyn_HTE_Lattice
    

    if geometry == "chain" #shortcut for chain
        lattice,LatGraph = get_finite_Lattice(2*L+1,"chain"; PBC = false)
        center_sites = [L+1]
        lattice.sitePositions = [lattice.sitePositions[i] .- lattice.sitePositions[L+1] for i in 1:2*L+1] #shift center site to zero coordinate
        return Dyn_HTE_Lattice(geometry ,lattice, LatGraph, center_sites)
    end

    ## helper function
    function replace_indices(tuple_vector,without)
        # Extract all unique integers from the tuples and sort them
        unique_values = unique(collect(Iterators.flatten(tuple_vector)))
        sorted_values = sort([x for x in unique_values if x ∉ without])
        
        # Create a mapping from each value to its ascending order index
        value_to_index = Dict(value => i for (i, value) in enumerate(sorted_values))
        
        # Replace each tuple entry with its corresponding index
        replaced_tuples = [tuple([if x ∉ without value_to_index[x] else i end for x in t ]...) for (i,t) in enumerate(tuple_vector)]
        
        return replaced_tuples
    end


    # Get the lattice and its corresponding graph representation
    lattice, LatGraph = get_finite_Lattice(2 * L + 1, geometry; PBC = false)

    # Extract the number of sites in the unit cell
    basis = length(lattice.unitcell.basis)

    # Compute the indices of the "center" vertices in the lattice
    center_vertices = [basis * sum([(2L + 1)^n * L for n = 0:(length(lattice.size) - 1)]) + b for b in 1:basis]

    # Get the position of the first center vertex (used as a reference point later)
    center_pos = lattice.sitePositions[center_vertices[1]]

    # Calculate the shortest path distances from each center vertex to all other vertices
    distances = [dijkstra_shortest_paths(LatGraph, center_vertices[i]).dists for i = 1:basis]

    # Identify vertices that are farther away than a threshold distance `L`
    too_large_b = [findall(x -> x > L, distances[i]) for i = 1:basis]

    # Find the common set of "too large" vertices across all basis distances
    too_large = reduce(intersect, [too_large_b[i] for i = 1:basis])

    # Remove the vertices that are too far from the lattice center
    deleteat!(lattice.sitePositions, too_large)
    lattice.length = length(lattice.sitePositions)  # Update the lattice length
    deleteat!(lattice.interactionSites, too_large)  # Remove interactions involving deleted sites

    # Re-center the remaining site positions relative to the initial center vertex
    lattice.sitePositions = [pos - center_pos for pos in lattice.sitePositions]

    # Update the indices of interaction sites after removing "too large" vertices
    lattice.interactionSites = replace_indices(lattice.interactionSites, too_large)

    # Convert the updated lattice back into a graph
    gg = latticeToGraph(lattice)

    # Find the central site(s) of the updated graph (e.g., for further analysis)
    center_sites = find_graph_center(gg)

    # Return the modified lattice, its graph representation, and the central site(s)
    return Dyn_HTE_Lattice(geometry ,lattice, gg, center_sites)
end


"""
    find_site_basis_label(lattice) -> Vector{Int}

Identify which basis site in the unit cell each lattice site corresponds to.

# Arguments
- `lattice`: Lattice object with unit cell information

# Returns
- `Vector{Int}`: Array where each element indicates the basis site index (1 to n_basis)
  for the corresponding lattice site

# Description
For lattices with multiple sites per unit cell (like honeycomb or kagome),
this function determines which of the basis sites each position in the lattice
corresponds to. It works by checking if a site position minus a basis site
position can be expressed as an integer combination of lattice vectors.

# Examples
```julia
# For a honeycomb lattice (2 basis sites)
labels = find_site_basis_label(honeycomb_lattice)
# Returns array with values 1 or 2 for each site
```
"""
function find_site_basis_label(lattice)
    function is_int_vec(vec)
        return all(isinteger,(x->round(x; digits = 10)).(vec))
        end
    
        site_basis_label = zeros(Int,lattice.length)
        basis = collect.(lattice.unitcell.basis)
        lattice_vectors = reduce(hcat,(collect.(lattice.unitcell.primitive)))

        for (site_index,site) in enumerate(collect.(lattice.sitePositions))
        for (i,b) in enumerate(basis)
            fractional_coords = lattice_vectors \ (site - b)
            if is_int_vec(fractional_coords)
                site_basis_label[site_index] = Int(i)
            end
        end
        end

    return site_basis_label
end



###### TESTS ############
#lattice,LatGraph = getLattice(4,"honeycomb"; PBC = false);
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

#lattice,LatGraph,center_sites = getLattice_Ball(6,"honeycomb");
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))


# Function to determine which basis site a given point corresponds to
