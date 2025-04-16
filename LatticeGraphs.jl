using Graphs, GraphRecipes, Plots
#include Lattice support from SpinMC.jl
include("Lattice.jl")
include("LatticeSymmetries.jl")

function latticeToGraph(lattice::Lattice)::SimpleGraph{Int}
    """ transforms a Lattice Struct into a SimpleGraph """
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

function get_finite_Lattice(L::Int,geometry::String; PBC::Bool = true)
    """ creates lattice and corresponding graphs, L is the linear size, PBC sets the use of boundary conditions.
    Currently implmented:
        - chain
        - square
        - triang
        - honeycomb
        - pyrochlore
        - kagome
    """
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

### helper functions
import Base.-
function -(a::Tuple{Float64, Float64},b::Tuple{Float64, Float64})::Tuple{Float64, Float64}
    """ - for 2-tuples """
    return (a[1]-b[1],a[2]-b[2])
end
function -(a::Tuple{Float64, Float64, Float64},b::Tuple{Float64, Float64, Float64})::Tuple{Float64, Float64, Float64}
    """ minus for 3-tuples """
    return (a[1]-b[1],a[2]-b[2],a[3]-b[3])
end

function find_graph_center(graph)
    """ 
    Gives the sites at the center of the graph
    """

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

function getLattice(L::Int,geometry::String)::Dyn_HTE_Lattice
    """ 
    Gives the lattice where all sites are at most L away from the center sites, no PBC
    """

    if geometry == "chain" #shortcut for chain
        lattice,LatGraph = get_finite_Lattice(2*L+1,"chain"; PBC = true)
        center_sites = [L+1]
        return lattice,LatGraph,center_sites
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

###### TESTS ############
#lattice,LatGraph = getLattice(4,"honeycomb"; PBC = false);
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

#lattice,LatGraph,center_sites = getLattice_Ball(6,"honeycomb");
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))


# Function to determine which basis site a given point corresponds to
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
