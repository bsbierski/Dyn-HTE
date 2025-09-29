### find embedding factors

"""
is_simple_isomorphic(gG1::GraphG, gG2::GraphG) -> Bool

Check if the underlying simple graphs of two GraphG objects are isomorphic.

# Arguments
- `gG1::GraphG`: First graph with external legs. See [`GraphG`](@ref).
- `gG2::GraphG`: Second graph with external legs

# Returns
- `Bool`: `true` if the underlying simple graphs are isomorphic, `false` otherwise

# Description
Determines whether the two graphs have the same topology when edge weights 
are ignored, while respecting the external legs. This function handles special 
cases where one graph has coincident external legs (j=j') and the other doesn't.

# Examples
```julia
if is_simple_isomorphic(graph1, graph2)
    println("These graphs have the same underlying structure")
end
```
"""
function is_simple_isomorphic(gG1::GraphG,gG2::GraphG)::Bool
    ### convert gg,gg_flip to SimpleGraphs
    gg1_simple = toSimpleGraph(gG1.g)

    ### catch the cases of one gG being onsite and the other not 
    if gG1.jjp[1] == gG1.jjp[2] && gG2.jjp[1] != gG2.jjp[2] 
        return 0
    end

    if gG1.jjp[1] != gG1.jjp[2] && gG2.jjp[1] == gG2.jjp[2] 
        return 0
    end

    gg1_simple_color_vec = zeros(Int64,nv(gG1.g))
    gg2_simple_color_vec = zeros(Int64,nv(gG2.g))

    if gG1.jjp[1] == gG1.jjp[2] 
        if gG2.jjp[1] == gG2.jjp[2] 
        gg1_simple_color_vec[gG1.jjp[1]] = 1
        gg2_simple_color_vec[gG2.jjp[1]] = 1
        else
            return false
        end
    else
        gg1_simple_color_vec[[gG1.jjp[1],gG1.jjp[2]]] = [1,1]
        gg2_simple_color_vec[[gG2.jjp[1],gG2.jjp[2]]] = [1,1]
    end

    vertex_relation(j,i) = (gg1_simple_color_vec[j] == gg2_simple_color_vec[i])

    return has_isomorph(gg1_simple,gG2.g ; vertex_relation = vertex_relation)
end

###Initialize the File if it does not exist yet.
#vector = [0,[[gG_vec[1][1],[[0,1,1,true]],0]]]
#@save data_dir()*"/GraphFiles/unique_gG_vec_0.jld2" vector

"""
    give_unique_gG_vec(gG_vec::Vector{Vector{GraphG}}) -> unique_Graphs

Process graphs to identify unique topological representatives and group equivalent graphs.

# Arguments
- `gG_vec::Vector{Vector{GraphG}}`: Nested vector of graphs with external legs organized by order

# Returns
- `unique_Graphs`: Structure containing unique graph representatives and their equivalent graphs.
  See [`unique_Graphs`](@ref) for the complete structure definition.

# Description
This function processes a collection of `GraphG` objects to identify those with equivalent underlying
simple graph topology (ignoring edge weights), grouping them together for efficient
computation. This optimization is crucial for embedding calculations, as graphs with
the same topology have identical embedding factors.

The function creates a `unique_Graphs` structure containing:
- `max_order::Int`: The highest order processed
- `graphs::Vector{unique_Graph}`: Array of unique graph representatives

Each `unique_Graph` contains:
- `ref_graph::GraphG`: A representative graph from the equivalence class
- `distance::Int`: Graph distance between external vertices
- `gG_vec::Vector{gG_properties}`: Properties of all equivalent graphs

The function implements caching:
1. First attempts to load precomputed results from disk
2. If not found, loads results for one order lower and extends incrementally
3. Progress is saved to disk to avoid recomputation

For each new graph, the algorithm:
1. Checks if it's topologically equivalent to any existing unique representative using [`is_simple_isomorphic`](@ref)
2. If equivalent, adds its `gG_properties` to that representative's group
3. If unique, creates a new `unique_Graph` entry with distance information

# Performance Notes
- Processes graphs in batches of 1000 with progress reporting
- Results are cached to `/GraphFiles/unique_gG_vec_<maxorder>.jld2`

# Examples
```julia
# Process all graphs up to order 10
unique_graphs = give_unique_gG_vec(all_graphs_by_order)
@assert unique_graphs isa unique_Graphs

# Use in embedding calculations
coeffs = Calculate_Correlator_fast(lattice, i, j, 10, unique_graphs, graph_coeffs)
```

# See Also
- [`unique_Graphs`](@ref): The returned data structure
- [`unique_Graph`](@ref): Individual equivalence classes of graphs
- [`gG_properties`](@ref): Properties stored for each equivalent graph
- [`is_simple_isomorphic`](@ref): Determines topological equivalence between graphs
- [`Calculate_Correlator_fast`](@ref): Uses unique graph structure for efficient correlator calculation
"""
function give_unique_gG_vec(gG_vec::Vector{Vector{GraphG}})

    
    maxorder = length(gG_vec) - 1

    # try to load the file. if it does not exist try to load the file of one less order
    file_path = data_dir()*"/GraphFiles/unique_gG_vec_$maxorder.jld2"

    if isfile(file_path)
        unique_gG_vec = load_object(file_path) 
        return unique_gG_vec
    else
        # if the order has not been calculated: try one order less
        unique_gG_vec = give_unique_gG_vec(gG_vec[1:end-1])
    end

    unique_order = unique_gG_vec[1]
    
    if maxorder <= unique_order
        if maxorder == unique_order
            return  unique_gG_vec
        else
            #todo delete all graphs of order larger than "maxorder"
            return  unique_gG_vec
        end
    else
        #unique_gGs = map(x -> x[1], unique_gG_vec[2])
        for (o,gG_vec_order) in enumerate(gG_vec[(unique_order+2):end])
            for (index,gg) in enumerate(gG_vec_order)
                if mod(index,1000) ==0
                    println(string(index)*" out of "*string(length(gG_vec_order))*" done")
                end
                # find first isomorphic graph to gg that is already in the unique list. There is at most one. 
                unique_index = findfirst(x->is_simple_isomorphic(gg,x[1]), unique_gG_vec[2])
                # if there is no matching graph: add gg to list of unique graphs
                if unique_index === nothing
                    dist = dijkstra_shortest_paths(gg.g,gg.jjp[1]).dists[gg.jjp[2]]
                    push!(unique_gG_vec[2],[ gg,[[unique_order + o  ,index,symmetryFactor(gg),is_symmetric(gg)]],dist])
                else
                    push!(unique_gG_vec[2][unique_index[1]][2],[unique_order + o  ,index,symmetryFactor(gg),is_symmetric(gg)])
                end

            end
        end
    end

    unique_gG_vec[1] = maxorder
    @save data_dir()*"/GraphFiles/unique_gG_vec_$maxorder.jld2" unique_gG_vec

    return unique_gG_vec

end

"""
    give_unique_gG_vec(max_order::Int) -> unique_Graphs

Load precomputed unique graph representatives up to the specified maximum order.

# Arguments
- `max_order::Int`: Maximum order of graphs to load

# Returns
- `unique_Graphs`: Structure containing unique graph representatives and their equivalent graphs.
  See [`unique_Graphs`](@ref) for the complete structure definition.

# Description
Loads precomputed `unique_Graphs` structures from disk for the specified maximum order.
This function provides access to cached results from previous computations, avoiding
the need to reprocess large collections of `GraphG` objects.

The returned `unique_Graphs` structure contains:
- `max_order::Int`: Maximum order of graphs in the collection  
- `graphs::Vector{unique_Graph}`: Array of unique graph representatives

Each `unique_Graph` in the collection provides:
- `ref_graph::GraphG`: Representative graph for the equivalence class
- `distance::Int`: Graph distance between external vertices
- `gG_vec::Vector{gG_properties}`: Properties of all graphs equivalent to the representative.

The function handles several scenarios:
1. **Direct load**: If `unique_gG_vec_<max_order>.jld2` exists, loads it directly
2. **Error handling**: Throws `ArgumentError` if no suitable file is found

# File Structure
The function looks for files in the pattern:
- `GraphFiles/unique_gG_vec_<order>.jld2` for most orders
- `GraphFiles/unique_gG_vec_12_<part>.jld2` for order 12 parts (contains `unique_graphs_12_part`)

# Throws
- `ArgumentError`: If no suitable graph file is available for the requested order

# Examples
```julia
# Load unique graph structure up to order 10
unique_graphs = give_unique_gG_vec(10)
@assert unique_graphs isa unique_Graphs
@assert unique_graphs.max_order == 10

# Handle the special case of order 12 (automatically combines parts)
unique_graphs_12 = give_unique_gG_vec(12)

# Access individual unique graphs
for unique_graph in unique_graphs.graphs
    ref = unique_graph.ref_graph        # GraphG representative
    dist = unique_graph.distance        # Int distance
    props = unique_graph.gG_vec        # Vector{gG_properties}
end

# Use in correlator calculations
coeffs = Calculate_Correlator_fast(lattice, i, j, max_order, unique_graphs, coeffs)
```

# See Also
- [`unique_Graphs`](@ref): The returned data structure
- [`unique_Graph`](@ref): Individual equivalence classes within the structure  
- [`gG_properties`](@ref): Properties stored for equivalent graphs
- [`give_unique_gG_vec(::Vector{Vector{GraphG}})`](@ref): Processes raw graphs to create unique structure
- [`Calculate_Correlator_fast`](@ref): Main consumer of unique graph structures
"""
function give_unique_gG_vec(max_order::Int)
    # try to load the file. if it does not exist try to load the file of one less order
    file_path = data_dir()*"/GraphFiles/unique_gG_vec_$max_order"*".jld2"

    if isfile(file_path)
        unique_gG_vec = load_object(file_path) 
        return unique_gG_vec
    else
        if max_order == 12
                ##Combine the 4 parts for order 12
                graphlist = Vector{Vector{unique_Graph}}(undef,4)
                for part = 1:4
                    @load data_dir()*"/GraphFiles/unique_gG_vec_$max_order"*"_$part"*".jld2" unique_graphs_12_part
                    graphlist[part] = unique_Graphs(12,unique_graphs_12_part).graphs
                end
                    combined = vcat(graphlist...) 
                    combined_unique = unique_Graphs(max_order,combined)
                    @save data_dir()*"/GraphFiles/unique_gG_vec_12.jld2" combined_unique
            return combined_unique
        end
        throw(ArgumentError("No unique graph file available for order $max_order"))
    end

end


"""
    e_fast(LL::SimpleGraph{Int}, j::Int, jp::Int, gG::GraphG) -> Int

Calculate the embedding factor of a graph with external legs into a lattice.

# Mathematical Background
The lattice coefficients in the high-temperature expansion are related to the embedding 
factors and graph coefficients through:
```math
c_{ii'}^{(n)}(i\\nu_m) = \\sum_{g^{(n)}} e(\\mathcal{L},i,i',g^{(n)}) \\cdot c_{g^{(n)}}(i\\nu_m)
```
where:
- ``c_{ii'}^{(n)}(i\\nu_m)`` are the lattice coefficients for sites i, i' at order n
- ``e(\\mathcal{L},i,i',g^{(n)})`` is the embedding factor (computed by this function)
- ``c_{g^{(n)}}(i\\nu_m)`` are the coefficients on individual graphs g^{(n)}

This function computes the embedding factor ``e(ℒ,i,i',g^{(n)})`` for a specific graph.

# Arguments
- `LL::SimpleGraph{Int}`: The lattice graph ℒ
- `j::Int`: First external site in the lattice (i)
- `jp::Int`: Second external site in the lattice (i')
- `gG::GraphG`: Graph g^{(n)} with external legs to embed

# Returns
- `Int`: Number of ways the graph can be embedded with j, jp as external sites

# Description
Counts the number of isomorphic subgraphs in the lattice that match the given
graph structure, with the constraint that the external legs must map to the
specified lattice sites j and jp. This embedding factor quantifies how many
ways a particular graph topology can be realized in the lattice with the
specified external connections.

# Note
This is an optimized version that directly counts subgraph isomorphisms.

# Examples
```julia
# Count embeddings of a correlation graph between sites 1 and 5
embeddings = e_fast(lattice_graph, 1, 5, correlation_graph)
```
"""
function e_fast(LL::SimpleGraph{Int},j::Int,jp::Int,gG::GraphG)::Int
    numSubIsos = count_subgraphisomorph(LL,gG.g,jL1 = j,jL2 = jp,jG1 = gG.jjp[1],jG2 = gG.jjp[2])
    return numSubIsos 
end


"""
    Calculate_Correlator_fast(L::SimpleGraph{Int}, ext_j1::Int, ext_j2::Int, 
                             max_order::Int, gG_vec_unique::unique_Graphs, 
                             C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}}) 
                             -> Vector{Vector{Rational{Int128}}}

Calculate the coefficients of ``(-x)^n`` for ``TG_{j_1 j_2}(iν_m)`` using optimized embedding calculations.

# Mathematical Background
The thermal Green's function is expanded in powers of x = J/T as:
```math
TG_{j_1 j_2}(i\\nu_m) = \\sum_{n=0}^\\infty (-x)^n c_{j_1 j_2}^{(n)}(i\\nu_m), \\quad (x = \\frac{J}{T})
```
This function computes the coefficients ``c_{j_1 j_2}^{(n)}(iν_m)`` for the high-temperature series expansion.

# Arguments
- `L::SimpleGraph{Int}`: The lattice graph
- `ext_j1::Int`: First external site in the lattice 
- `ext_j2::Int`: Second external site in the lattice 
- `max_order::Int`: Maximum expansion order n to calculate
- `gG_vec_unique::unique_Graphs`: Structure of unique graph representatives
- `C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}}`: Coefficients for each graph

# Returns
- `Vector{Vector{Rational{Int128}}}`: Array of coefficients for each order and power of  ``\\Delta = \\frac{1}{\\nu_m}``

# Description
Computes the coefficients of the high-temperature series expansion for a correlation
function between sites ext_j1 and ext_j2. Each order returns a vector representing 
prefactors of [δω, Δ², Δ⁴, Δ⁶, Δ⁸, Δ¹⁰, Δ¹², Δ¹⁴, Δ¹⁶, Δ¹⁸].

The computation uses preprocessing of [`unique_Graphs`](@ref) structures to minimize redundant 
calculations of embedding factors. This optimization groups topologically equivalent graphs
into [`unique_Graph`](@ref) equivalence classes, allowing embedding factors to be computed
only once per unique topology via [`e_fast`](@ref).

# Examples
```julia
 # Calculate correlation coefficients between sites 1 and 5 up to order 10
 coeffs = Calculate_Correlator_fast(lattice_graph, 1, 5, 10, unique_graphs, graph_coeffs)
```
"""
function Calculate_Correlator_fast(L::SimpleGraph{Int},ext_j1::Int,ext_j2::Int,max_order::Int,gG_vec_unique::unique_Graphs,C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}})::Vector{Vector{Rational{Int128}}}

    #initialize result array
    result_array = Vector{Vector{Rational{Int128}}}(undef, max_order+1)

    #for every order we get result vector representing prefactors of [δw,Δ^2,Δ^4,Δ^6,Δ^8,Δ^10,Δ^12,Δ^14,Δ^16,Δ^18]
    for ord = 1:max_order+1
        result_array[ord] = zeros(Rational{Int128},10)
    end

    #calculate the shortest graph distance between ext_j1 and ext_j2
    ext_dist = dijkstra_shortest_paths(L,ext_j1).dists[ext_j2]

    # only iterate over the unique simple graphs in unique_Gg
    for (index,unique_Gg) in enumerate(gG_vec_unique.graphs)
        gg = unique_Gg.ref_graph   #Graph
        gg_dist = unique_Gg.distance #edge distance between the external vertices
          # if the graph is long enough
          if gg_dist < ext_dist 
            continue
        end

        #if its the onsite correlator we only need on-site graphs 
        if ext_dist == 0
            if gg_dist > ext_dist 
                continue
            end
        else #if not, we dont need any on site graphs
            if gg_dist == 0 
                continue
            end
        end

        #calculate the embedding factor
        
        

        emb_fac = e_fast(L,ext_j1,ext_j2,gg)

        emb_fac_assym = 0
        if all(x -> !x.is_symmetric, unique_Gg.gG_vec)
            emb_fac_assym = e_fast(L,ext_j2,ext_j1,gg)
        else
            emb_fac_assym = emb_fac
        end

       # println("$index th graph embeding factor = $emb_fac")

        #### now we sum overall graphG eqivalent to the unique Gg
        for graph in unique_Gg.gG_vec
            g_order = graph.order #order
            gG_vec_index = graph.index #index
            symmetry_factor = graph.symmetry_factor#symmetry factor
            is_symmetric = graph.is_symmetric  #bool if the graph is symmetric
           
            fac = emb_fac + emb_fac_assym
            if  is_symmetric
                fac = emb_fac
            end

            #look up the value of the graph from C_Dict_vec
            look_up_dict = C_Dict_vec[g_order+1][gG_vec_index]
            
            
            result_array[g_order+1] .+= look_up_dict.*Int128(fac/symmetry_factor)
        end
    end

    return result_array
end

###### LEGACY FUNCTIONS (SLOW)
"""
OLD FUNKTION: use e_fast
find embedding factor e of GraphG gG in lattice L
- lattice L (needs to be chosen large enough to avoid boundary effects!)
- external site indices jjp=[j,j′] can be (i,i′) or (i′,i) if gG is not symmetric under exchange of i ↔ i′
- assumes that the distance j-jp is smaller or equal to the distance of external vertices of gG
"""
function e(L::SimpleGraph{Int},j::Int,jp::Int,gG::GraphG)::Int

 
    ### copy L -> LL and gG.g -> gg not to mess with input
    LL = copy(L)
    gg = copy(gG.g)
    gg_flip = copy(gG.g)

    ### add the terminals j,j' to LL as bonds to extra vertices
    add_vertex!(LL)
    add_edge!(LL,j ,nv(LL))
    add_vertex!(LL)
    add_edge!(LL,jp,nv(LL))

    ### add to gg terminals i,i' and i',i with bond-weight 100,200
    add_vertex!(gg)
    add_edge!(gg,gG.jjp[1],nv(gg),100)
    add_vertex!(gg)
    add_edge!(gg,gG.jjp[2],nv(gg),200)

    add_vertex!(gg_flip)
    add_edge!(gg_flip,gG.jjp[2],nv(gg_flip),100)
    add_vertex!(gg_flip)
    add_edge!(gg_flip,gG.jjp[1],nv(gg_flip),200)

    ### convert gg,gg_flip to SimpleGraphs
    gg_simple = toSimpleGraph(gg)
    gg_flip_simple = toSimpleGraph(gg_flip)

    ### define rules for mapping i,i' to j,j'
    LL_color_vec = zeros(nv(LL))
    gg_simple_color_vec = zeros(nv(gg_simple))

    ### prepare isomorphism check on gg <--> ggflip respecting the edge-weights for the next if-clause
    edge_relation(b1,b2) = (gg.weights[src(b1),dst(b1)] == gg_flip.weights[src(b2),dst(b2)])

    ### if j=j' (-> i=i') or the graph gG,i,i' is isomorph to gG,i',i (respecting edge-weights)
    if j==jp || Graphs.Experimental.has_isomorph(gg_simple,gg_flip_simple,edge_relation=edge_relation)
        LL_color_vec[end-1:end] = [1,2]
        gg_simple_color_vec[end-1:end] = [1,2]
    
    ### otherwise need to allow both i,i' <--> j,j' and i,i' <--> j',j 
    else
        LL_color_vec[end-1:end] = [1,1]
        gg_simple_color_vec[end-1:end] = [1,1]
    end

    vertex_relation(j,i) = (LL_color_vec[j] == gg_simple_color_vec[i])

    numSubIsos = Graphs.Experimental.count_subgraphisomorph(LL,gg_simple,vertex_relation=vertex_relation)
    return numSubIsos / symmetryFactor(gG)
end


 """OLD FUNCTION: use Calculate_Correlator_fast"""
function Calculate_Correlator(L::SimpleGraph{Int},ext_j1::Int,ext_j2::Int,max_order,gG_vec::Vector{Vector{GraphG}},C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}})::Vector{Vector{Rational{Int128}}}
   

    #initialize result array
    result_array = Vector{Vector{Rational{Int128}}}(undef, max_order+1)


    #for every order we get result vector representing Delta^2 prefactors
    for ord = 1:max_order+1
        result_array[ord] = zeros(Rational{Int128},10)
    end


    #now calculate order for order the correlator

    Threads.@threads for gG_arr in gG_vec[1:max_order+1]

        for gG_idx in eachindex(gG_arr)
            g_order = Int(sum(gG_arr[gG_idx].g.weights)/2)

            #now we sum over all graphG

            look_up_dict = C_Dict_vec[g_order+1][gG_idx]

            emb_fac = e(L,ext_j1,ext_j2,gG_arr[gG_idx])


            result_array[g_order+1] += look_up_dict*emb_fac   
        end
    end


    return result_array
end
