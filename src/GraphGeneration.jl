### all tools revolving around abstract graphs and their visualization

plt_empty = plot(label="",axis=([], false))
graphsInRow = 6 #for plotting

####### various helper functions ################### 
"""
    load_dyn_hte_graphs(spin_length::Number, max_order::Int; verbose = false) -> Dyn_HTE_Graphs

Load pre-computed graph data for high-temperature series expansion calculations.

# Arguments
- `spin_length::Number`: The spin length (e.g., 1/2 or 1)
- `max_order::Int`: Maximum expansion order to load
- `verbose::Bool=false`: Whether to print loading information

# Returns
- `Dyn_HTE_Graphs`: Structure containing unique graphs and their coefficients

# Description
Loads the required graph data files and evaluation coefficients needed for 
high-temperature series expansion calculations. If files are split across multiple 
parts, they are automatically merged.

# Examples
```julia
# Load graph data for spin-1/2 system up to 10th order
graphs_data = load_dyn_hte_graphs(1/2, 10)
```
"""
function load_dyn_hte_graphs(spin_length::Number,max_order::Int; verbose = false)::Dyn_HTE_Graphs
    #first check if all required data files are available, if not merge them
    merge_data_files()

    S = rationalize(spin_length)
    if S == 1//2
        S_string = "Spin_S1half"
    elseif S == 1//1
        S_string = "Spin_S1"
    else
        throw(error("Spinlength "*string(spin_length)*" is not yet implemented"))
    end

    #load list of unique graphs
    gG_vec_unique = give_unique_gG_vec(max_order);

    #create vector of all lower order dictionaries
    C_Dict_vec = Vector{Vector{Vector{Rational{Int128}}}}(undef,max_order+1) ;
    #load dictionaries of all lower orders C_Dict_vec 
    for ord = 0:max_order
        if verbose
        println("loading C_"*string(ord)*" for S="*S_string)
        end
        C_Dict_vec[ord+1]  = load_object(data_dir()*"/GraphEvaluations/"*S_string*"/C_"*string(ord)*".jld2")
    end 
    Dyn_HTE_Graphs(S,gG_vec_unique,C_Dict_vec)
end

"""
    split_vec(vec::Vector, part::Int, parts::Int) -> (Vector, Int, Int)

Split a vector into parts and return a specific chunk along with its indices.

# Arguments
- `vec::Vector`: The vector to split
- `part::Int`: The part number to return (1-indexed)
- `parts::Int`: Total number of parts to split into

# Returns
- `Tuple{Vector, Int, Int}`: The requested chunk and its start/end indices

# Description
Divides a vector into roughly equal parts, with any remainder added to the last part.
Returns the requested part along with its starting and ending indices in the original vector.

# Examples
```julia
data = collect(1:100)
chunk, start_idx, end_idx = split_vec(data, 2, 4)  # Get second quarter of data
```
"""
function split_vec(vec::Vector,part::Int,parts::Int)
    
    chunkLen = Int(floor(length(vec)/parts))

    if parts==1 && part==1
        return vec,1,length(vec)
    end

    @assert part>0
    @assert parts>=part
    if part==parts
        return vec[1+(part-1)*chunkLen:end], 1+(part-1)*chunkLen, length(vec)
    else
        return vec[1+(part-1)*chunkLen:(part)*chunkLen], 1+(part-1)*chunkLen, (part)*chunkLen
    end
end


"""
    degeneracy(g::SimpleWeightedGraph{Int64, Int64}) -> Int

Calculate the degeneracy factor of a weighted graph.

# Arguments
- `g::SimpleWeightedGraph{Int64, Int64}`: The weighted graph

# Returns
- `Int`: Degeneracy factor as product of factorials of edge multiplicities

# Description
Computes the product of the factorials of edge multiplicities for all edges in the graph.
This factor is used in statistical weight calculations for graph-based expansions.

# Examples
```julia
g = SimpleWeightedGraph(3)
add_edge!(g, 1, 2, 2)  # Edge with weight 2
add_edge!(g, 2, 3, 1)  # Edge with weight 1
deg = degeneracy(g)    # Returns 2! * 1! = 2
```
"""
function degeneracy(g::SimpleWeightedGraph{Int64, Int64})::Int
    
    mat = g.weights
    deg = 1
    n = nv(g)
    for row in 1:n
        for col in row+1:n
            deg *= factorial(mat[row,col])
        end
    end
    return deg
end

"""
    isIsomorph(g1::Graph, g2::Graph) -> Bool

Check if two Graph objects are isomorphic, respecting edge weights.

# Arguments
- `g1::Graph`: First graph
- `g2::Graph`: Second graph

# Returns
- `Bool`: `true` if graphs are isomorphic, `false` otherwise

# Description
Tests whether there exists a vertex mapping between g1 and g2 that preserves
both the graph structure and the edge weights.

# Examples
```julia
if isIsomorph(graph1, graph2)
    println("These graphs are structurally equivalent")
end
```
"""
function isIsomorph(g1::Graph,g2::Graph)::Bool
   
    ### convert g1,2 to SimpleGraphs
    g1_simple = toSimpleGraph(g1.g)
    g2_simple = toSimpleGraph(g2.g)

    ### prepare isomomorphism check on g1,2_simple respecting edge weigths in g1,g2 
    edge_relation(b1,b2) = (g1.g.weights[src(b1),dst(b1)] == g2.g.weights[src(b2),dst(b2)])
    
    return Graphs.Experimental.has_isomorph(g1_simple,g2_simple,edge_relation=edge_relation)
end

"""
    isIsomorph(gG1::GraphG, gG2::GraphG) -> Bool

Check if two GraphG objects (graphs with external legs) are isomorphic.

# Arguments
- `gG1::GraphG`: First graph with external legs
- `gG2::GraphG`: Second graph with external legs

# Returns
- `Bool`: `true` if graphs are isomorphic, `false` otherwise

# Description
Tests whether there exists a vertex mapping between gG1 and gG2 that preserves
the graph structure, edge weights, and the relationship between external legs.
External legs can match either in the same order or flipped.

# Examples
```julia
if isIsomorph(graph_with_legs1, graph_with_legs2)
    println("These graphs with external legs are equivalent")
end
```
"""
function isIsomorph(gG1::GraphG,gG2::GraphG)::Bool
    
    
    ### add terminal vertices via edge with weight 100 
    gG1_ext = copy(gG1.g)
    gG2_ext = copy(gG2.g)
    
    for j in 1:2
        add_vertex!(gG1_ext)
        add_edge!(gG1_ext,gG1.jjp[j],nv(gG1_ext),100)

        add_vertex!(gG2_ext)
        add_edge!(gG2_ext,gG2.jjp[j],nv(gG2_ext),100)
    end
    
    ### convert g1G,gG2 to SimpleGraphs
    gG1_ext_simple = toSimpleGraph(gG1_ext)
    gG2_ext_simple = toSimpleGraph(gG2_ext)

    ### prepare isomomorphism check on g1,2_simple respecting edge weigths in g1,g2 
    edge_relation(b1,b2) = (gG1_ext.weights[src(b1),dst(b1)] == gG2_ext.weights[src(b2),dst(b2)])
    
    return Graphs.Experimental.has_isomorph(gG1_ext_simple,gG2_ext_simple,edge_relation=edge_relation)
end

"""
find if gG is a symmetric graph with respect to switching the two external legs (fix the mapping of the two external vertices)
"""
function is_symmetric(gG::GraphG)::Bool

    gg = gG.g
    gg_simple = toSimpleGraph(gg)

    edge_relation(b1,b2) = (gg.weights[src(b1),dst(b1)] == gg.weights[src(b2),dst(b2)])

    # finds if there is an isomorphism by only permuting the internal vertices between the graph and the graph with its external vertices flipped.
    count = count_subgraphisomorph(gg_simple,gg,edge_relation=edge_relation,jL1 = gG.jjp[2],jL2 = gG.jjp[1],jG1 = gG.jjp[1],jG2 = gG.jjp[2])
    
    if count >0 
        return true
    else
        return false
    end
end

 """ search g in g_vec, if found return index, else return 0"""
function findg(g::Graph,g_vec::Vector{Graph})::Int
   
    for k in eachindex(g_vec)
        if isIsomorph(g,g_vec[k])
            return k
        end
    end
    return 0
end

""" find gG in gG_vec, if found return index, else return 0"""
function findg(gG::GraphG,gG_vec::Vector{GraphG})::Int
    
    for k in eachindex(gG_vec)
        if isIsomorph(gG,gG_vec[k])
            return k
        end
    end
    return 0
end

""" take Graph g and add one edge in all possible ways (strenghten existing edge or connect to new vertex) """
function addOneEdge(g::Graph)::Vector{Graph}
    
    vs = nv(g.g)
    g_vec = Graph[]

    ### add edge between any existing vertices 
    for src in 1:vs-1, dest in src+1:vs
        gnew = copy(g.g)
        add_edge!(gnew,src,dest,get_weight(g.g,src,dest)+1)
        push!(g_vec,Graph(gnew))
    end

    ### or create new vertex with edge to any existing vertex
    for src in 1:vs
        gnew = copy(g.g)
        add_vertex!(gnew)
        add_edge!(gnew,src,vs+1)
        push!(g_vec,Graph(gnew))
    end

    return g_vec
end

"""
    totalEdges(g::SimpleWeightedGraph{Int64, Int64}) -> Int

Count the total number of edges in a weighted graph.

# Arguments
- `g::SimpleWeightedGraph{Int64, Int64}`: The weighted graph

# Returns
- `Int`: The total number of edges, counting multiplicities

# Description
Computes the sum of all edge weights divided by 2 to account for the 
undirected nature of the graph (each edge is counted twice in the adjacency matrix).

# Examples
```julia
g = SimpleWeightedGraph(3)
add_edge!(g, 1, 2, 2)  # Double edge between 1 and 2
add_edge!(g, 2, 3, 1)  # Single edge between 2 and 3
total = totalEdges(g)  # Returns 3
```
"""
function totalEdges(g::SimpleWeightedGraph{Int64, Int64})::Int
    return Int(sum(g.weights)/2)
end


""" do not affect trivial graph with single vertex """
function removeVerticesWithoutEdge!(g::Graph)
    
    for j in 2:nv(g)
        if length(outneighbors(g,j))==0
            rem_vertex!(g,j)
        end
    end
end

"""
    toSimpleGraph(g::SimpleWeightedGraph{Int64, Int64}) -> SimpleGraph

Convert a weighted graph to a simple unweighted graph.

# Arguments
- `g::SimpleWeightedGraph{Int64, Int64}`: The weighted graph to convert

# Returns
- `SimpleGraph`: An unweighted graph with the same connectivity

# Description
Creates a simple graph where edges exist if and only if the corresponding
edges in the weighted graph have positive weights. The resulting graph
preserves connectivity but discards weight information.

# Examples
```julia
g_weighted = SimpleWeightedGraph(3)
add_edge!(g_weighted, 1, 2, 2)  # Edge with weight 2
add_edge!(g_weighted, 2, 3, 1)  # Edge with weight 1
g_simple = toSimpleGraph(g_weighted)  # Simple graph with edges 1-2 and 2-3
```
"""
function toSimpleGraph(g::SimpleWeightedGraph{Int64, Int64})::SimpleGraph
    ### drop the weights of a SimpleWeightedGraph to get a SimpleGraph
    g_simple=SimpleGraph(nv(g))
    for v1 in 1:nv(g), v2 in v1+1:nv(g)
        if g.weights[v1,v2]>0
            add_edge!(g_simple,v1,v2)
        end
    end
    return g_simple
end


"""
    splitToConnectedComp(g::SimpleWeightedGraph{Int64, Int64}) -> Vector{Graph}

Split a graph into its connected components, discarding isolated vertices.

# Arguments
- `g::SimpleWeightedGraph{Int64, Int64}`: The graph to split

# Returns
- `Vector{Graph}`: Array of connected components as separate Graph objects

# Description
Decomposes a graph into its connected components, returning each component
as a separate Graph object. Isolated vertices (single-vertex components) are
discarded from the result.

# Examples
```julia
# Split a disconnected graph into its components
components = splitToConnectedComp(disconnected_graph)
```
"""
function splitToConnectedComp(g::SimpleWeightedGraph{Int64, Int64})::Vector{Graph}
   
    if is_connected(g) 
        return [Graph(g)]
    else
        res_vec = Graph[]
        for vs in connected_components(g)
            if length(vs)>1
                g_SWG_split = SimpleWeightedGraph(g.weights[vs,vs])
                push!(res_vec,Graph(g_SWG_split))
            end
        end
        return res_vec
    end
end

"""
    numberOfLeaves(g::SimpleWeightedGraph{Int64, Int64}) -> Int

Count the number of leaf vertices in a graph.

# Arguments
- `g::SimpleWeightedGraph{Int64, Int64}`: The graph to analyze

# Returns
- `Int`: Number of leaf vertices

# Description
A leaf vertex is defined as a vertex with total edge weight less than 2.
This function counts all such vertices in the graph.

# Examples
```julia
leaves = numberOfLeaves(my_graph)
```
"""
function numberOfLeaves(g::SimpleWeightedGraph{Int64, Int64})::Int
    ### get number of leaves
    res = 0
    for i in vertices(g)
        if sum(g.weights[i,:])<2
            res = res+1
        end
    end
    return res
end

"""
    noLeavesExceptAt(g::SimpleWeightedGraph{Int64, Int64}, j_vec::Vector{Int64}=Int64[]) -> Bool

Check if a graph has leaves only at specified vertices.

# Arguments
- `g::SimpleWeightedGraph{Int64, Int64}`: The graph to check
- `j_vec::Vector{Int64}=Int64[]`: Vertices where leaves are allowed

# Returns
- `Bool`: `true` if the graph has no leaves except possibly at vertices in `j_vec`

# Description
A leaf vertex is defined as a vertex with total edge weight less than 2.
This function checks that such vertices occur only in the specified list `j_vec`.

# Examples
```julia
# Check if graph has leaves only at vertices 1 and 5
if noLeavesExceptAt(g, [1, 5])
    println("Graph has no unexpected leaf vertices")
end
```
"""
function noLeavesExceptAt(g::SimpleWeightedGraph{Int64, Int64},j_vec::Vector{Int64}=Int64[])::Bool
    ### check if a Graph g has leaves, exclude sites in j_vec from checking
    for i in vertices(g)
        bonds_i = sum(g.weights[i,:])
        if bonds_i<2 && !(i in j_vec)
            return false
        end
    end
    return true
end

"""
    hasGeneralizedLeaves(g::Graph) -> Bool

Check if a graph has generalized leaves.

# Arguments
- `g::Graph`: The graph to check

# Returns
- `Bool`: `true` if the graph has generalized leaves, `false` otherwise

# Description
A generalized leaf is a bridge edge (weight 1) whose removal would disconnect
the graph. This function checks for the existence of such edges.

# Examples
```julia
if hasGeneralizedLeaves(my_graph)
    println("Graph has generalized leaves (bridge edges)")
end
```
"""
function hasGeneralizedLeaves(g::Graph)::Bool
    
    g_ext = copy(g.g)

    ### move through all edges e, if weigth=1 remove it and check if there is still just one connected component 
    for e in edges(g_ext)
        if get_weight(g_ext,src(e),dst(e))==1
            g_ext_rem = copy(g_ext)
            rem_edge!(g_ext_rem,e)
            g_ext_rem = toSimpleGraph(g_ext_rem)
            con_comp = connected_components(g_ext_rem)
            if length(con_comp)>1 
                return true
            end
        end
    end
    return false
end 


"""
    hasGeneralizedLeaves(gG::GraphG) -> Bool

Check if a graph with external legs has generalized leaves.

# Arguments
- `gG::GraphG`: The graph with external legs to check

# Returns
- `Bool`: `true` if the graph has generalized leaves, `false` otherwise

# Description
For graphs with external legs, a generalized leaf is a bridge edge (weight 1)
whose removal would create a separate component that still connects the
external vertices. This function checks for the existence of such edges.

# Examples
```julia
if hasGeneralizedLeaves(my_graph_with_legs)
    println("Graph has generalized leaves (bridge edges)")
end
```
"""
function hasGeneralizedLeaves(gG::GraphG)::Bool
    
    
    ### add terminal vertices via edge with weight 100 
    gG_ext = copy(gG.g)
    
    for j in 1:2
        add_vertex!(gG_ext)
        add_edge!(gG_ext,gG.jjp[j],nv(gG_ext),100)
    end

    ### move through all edges e, if weigth=1 remove it, if j,jp are still connected see if we get a second disconnected component 
    for e in edges(gG_ext)
        if get_weight(gG_ext,src(e),dst(e))==1
            gG_ext_rem = copy(gG_ext)
            rem_edge!(gG_ext_rem,e)
            gG_ext_rem = toSimpleGraph(gG_ext_rem)
            con_comp = connected_components(gG_ext_rem)
            if length(con_comp)>1 
                Nv = nv(gG.g)
                for vs in con_comp
                    if issubset([Nv+1,Nv+2],vs)
                        return true
                    end
                end
            end
        end
    end
    return false
end 

"""
    symmetryFactor(g::Graph) -> Int

Calculate the symmetry factor of a graph.

# Arguments
- `g::Graph`: The graph to analyze

# Returns
- `Int`: The symmetry factor (number of automorphisms)

# Description
Computes the number of graph automorphisms (self-isomorphisms) that preserve
edge weights. This factor is important for statistical weights in diagrammatic
expansions.

# Examples
```julia
sym = symmetryFactor(my_graph)
```
"""
function symmetryFactor(g::Graph)::Int
   
    gg = copy(g.g)

    ### convert gg to SimpleGraphs
    gg_simple = toSimpleGraph(gg)

    ### prepare auto-momorphism check on gg_simple respecting edge weigths in gg 
    edge_relation(b1,b2) = (gg.weights[src(b1),dst(b1)] == gg.weights[src(b2),dst(b2)])

    return Graphs.Experimental.count_isomorph(gg_simple,gg_simple,edge_relation=edge_relation)
end

"""
    symmetryFactor(gG::GraphG) -> Int

Calculate the symmetry factor of a graph with external legs.

# Arguments
- `gG::GraphG`: The graph with external legs to analyze

# Returns
- `Int`: The symmetry factor (number of automorphisms)

# Description
Computes the number of graph automorphisms (self-isomorphisms) that preserve
edge weights and the identity of external legs. This factor is important for
statistical weights in diagrammatic expansions with external sources.

# Examples
```julia
sym = symmetryFactor(my_graph_with_legs)
```
"""
function symmetryFactor(gG::GraphG)::Int
    
    gg = copy(gG.g)

    ### add to gg two vertices at terminals j,j' with bond-weight 100,101
    add_vertex!(gg)
    add_edge!(gg,gG.jjp[1],nv(gg),100)
    add_vertex!(gg)
    add_edge!(gg,gG.jjp[2],nv(gg),101)

    ### convert gg,gg_flip to SimpleGraphs
    gg_simple = toSimpleGraph(gg)

    ### prepare auto-momorphism check on gg_simple respecting edge weigths in gg 
    edge_relation(b1,b2) = (gg.weights[src(b1),dst(b1)] == gg.weights[src(b2),dst(b2)])

    return Graphs.Experimental.count_isomorph(gg_simple,gg_simple,edge_relation=edge_relation)
end

"""
    gplot(g::Graph; title::String="", save::Bool=false) -> Plot

Create a visualization of a Graph.

# Arguments
- `g::Graph`: The graph to plot
- `title::String=""`: Optional title for the plot
- `save::Bool=false`: Whether to save the plot to file

# Returns
- A plot object representing the graph visualization

# Description
Generates a visual representation of the graph with numbered vertices.
Edge weights are represented by multiple parallel edges.

# Examples
```julia
# Create and display a plot of a graph
plot = gplot(my_graph, title="Example Graph")
display(plot)

# Save the plot to a file
gplot(my_graph, title="saved_graph", save=true)
```
"""
function gplot(g::Graph;title::String="",save::Bool=false)
    Random.seed!(2011);
    m = collect(g.g.weights)
    n = size(m)[1]

    ### catch case of no bonds --- does not work yet!!!
    if maximum(m)==0
        fig=graphplot(g,names=1:n,markersize=0.22,fontsize=12,nodeshape=:circle)
    else
        ### need to convert in multigraph language that yields the destination sites
        gg =Vector{Int}[]
        for src in 1:n-1
            dests = Int[]
            for dest in src+1:n
                if m[src,dest]>0
                    append!(dests,[dest for _ in 1:m[src,dest]])
                end
            end
            append!(gg,[dests])
        end
        fig=graphplot(gg,names=[string(n) for n in 1:n],markercolor=:lightblue,markersize=0.15,fontsize=14,nodeshape=:circle,arrow=false,edgewidth=(s,d,w)->3,title=title)
    end

    ### save fig as png if filename supplied
    if save
        savefig("figs/"*title*".png")
    end
    return fig
end


"""
    gplot(gG::GraphG; title::String="", save::Bool=false) -> Plot

Create a visualization of a GraphG (graph with external legs).

# Arguments
- `gG::GraphG`: The graph with external legs to plot
- `title::String=""`: Optional title for the plot
- `save::Bool=false`: Whether to save the plot to file

# Returns
- A plot object representing the graph visualization

# Description
Generates a visual representation of the graph where external legs
are highlighted in red and labeled as "j" and "j'".

# Examples
```julia
# Create and display a plot of a graph with external legs
plot = gplot(my_graph_with_legs, title="Correlation Function Graph")
display(plot)
```
"""
function gplot(gG::GraphG;title::String="",save::Bool=false)
    Random.seed!(2011);
    gplot = copy(gG.g) #the SimpleWeightedGraph for plotting with external legs

    ### add external bonds in red
    for j in gG.jjp
        add_vertex!(gplot)
        add_edge!(gplot,j,nv(gplot))
    end

    m = collect(gplot.weights)
    n = size(m)[1]

    function edgewidth_function(s,d,w) 
        if d>=n-1
            return 1
        end
        return 3
    end

    ### need to convert in multigraph language that yields the destination sites
    gg =Vector{Int}[]
    for src in 1:n-1
        dests = Int[]
        for dest in src+1:n
            if m[src,dest]>0
                append!(dests,[dest for _ in 1:m[src,dest]])
            end
        end
        append!(gg,[dests])
    end
    fig=graphplot(gg,names=append!([string(n) for n in 1:n-2],["j","j'"]),markercolor=append!([:lightblue for _ in 1:n-2],[:red,:red]),nodesize=0.15,nodeweights=append!(ones(n-2),[0.4,0.4]),fontsize=14,nodeshape=append!([:circle for _ in 1:n-2],[:rect,:rect]),ew=edgewidth_function,arrow=false,title=title)

    ### save fig as png if filename supplied
    if save
        savefig("figs/"*title*".png")
    end
    return fig
end
"""
    gplot(g_vec::Vector; subtitle_vec::Vector{String}=["#"*string(pos) for pos in eachindex(g_vec)], 
          title::String="", save::Bool=false) -> Plot

Create a grid visualization of multiple graphs.

# Arguments
- `g_vec::Vector`: Collection of graphs to plot
- `subtitle_vec::Vector{String}`: Optional subtitles for each graph
- `title::String=""`: Optional title for the entire plot
- `save::Bool=false`: Whether to save the plot to file

# Returns
- A plot object representing the grid of graph visualizations

# Description
Generates a grid layout containing visualizations of multiple graphs,
with customizable subtitles for each graph and a main title.

# Examples
```julia
# Plot a collection of graphs in a grid
plot = gplot(graph_collection, title="Evolution of Graphs")
display(plot)

# Plot with custom subtitles
plot = gplot(graph_collection, subtitle_vec=["Initial", "Step 1", "Step 2", "Final"])
```
"""
function gplot(g_vec::Vector{};subtitle_vec::Vector{String}=["#"*string(pos) for pos in eachindex(g_vec)],title::String="",save::Bool=false)
    fig_vec = []
    for g_pos in 1:length(g_vec)
        fig = gplot(g_vec[g_pos];title=subtitle_vec[g_pos])
        append!(fig_vec,[fig])
    end
    rows = Int(floor(length(fig_vec)/graphsInRow))
    if rows*graphsInRow < length(g_vec) rows += 1 end
    fig = plot(fig_vec..., layout=(rows,graphsInRow),dpi=300,size=(300*graphsInRow,300*rows),plot_title=title)
    if save
        savefig("figs/"*title*".png")
    end
    return fig
end

"""
    getAllGraphsNextOrder(g_vec::Vector{Graph}) -> Vector{Graph}

Generate all connected graphs with one more edge than those in the input vector.

# Arguments
- `g_vec::Vector{Graph}`: Vector of graphs of order n

# Returns
- `Vector{Graph}`: All unique connected graphs of order n+1

# Description
Starting from a vector of graphs with n edges, generates all possible connected
graphs with n+1 edges, filtering out isomorphic duplicates. The resulting graphs
are limited to those with at most 2 leaf vertices.

# Examples
```julia
# Generate all order-5 graphs from order-4 graphs
graphs_order5 = getAllGraphsNextOrder(graphs_order4)
```

# Note
The resulting graphs are also saved to disk in the data directory.
"""
function getAllGraphsNextOrder(g_vec::Vector{Graph})::Vector{Graph}
    
    
    println("finding all graphs of order "*string(totalEdges(g_vec[end].g)+1)*"..." )
    
    g_new_vec = Graph[]
    len_g_vec = length(g_vec)
    for g_pos in eachindex(g_vec)
        g = g_vec[g_pos]
        println("adding on graph r="*string(g_pos)*" out of "*string(len_g_vec))
        g_test_vec = addOneEdge(g) #these are all unique
        for g_test in g_test_vec
            if numberOfLeaves(g_test.g) <= 2 ### && !is_cyclic(g_test.g) && maximum([length(neighbors(g_test.g, v)) for v in 1:nv(g_test.g)])<3 ### changes for chain_only graphs
            isNew = true
                Threads.@threads for g_new in g_new_vec   
                    if isIsomorph(g_new,g_test) 
                        isNew = false
                        break
                    end
                end 
                if isNew append!(g_new_vec,[g_test]) end
            end
        end
    end
    save_object(data_dir()*"/GraphFiles/graphs_"*string(totalEdges(g_new_vec[end].g))*".jld2",g_new_vec)
    return g_new_vec
end

"""
    getVacGraphs(graphs_vec::Vector{Vector{Graph}}) -> Vector{Vector{Graph}}

Filter a collection of graphs to remove those with generalized leaves.

# Arguments
- `graphs_vec::Vector{Vector{Graph}}`: Nested vector of graphs by order

# Returns
- `Vector{Vector{Graph}}`: Filtered collection with generalized leaves removed

# Description
From a collection of graphs organized by order, returns a new collection
where graphs with generalized leaves have been removed from all orders
except the first. This is used to identify "vacuum" graphs in diagrammatic
expansions.

# Examples
```julia
# Get vacuum graphs from a collection of all graphs
vac_graphs = getVacGraphs(all_graphs_by_order)
```
"""
function getVacGraphs(graphs_vec::Vector{Vector{Graph}})::Vector{Vector{Graph}}
    graphs_vac_vec = copy(graphs_vec)
    for m in 1:length(graphs_vac_vec)-1
        graphs_vac_vec[m+1] = [g for g in graphs_vac_vec[m+1] if !hasGeneralizedLeaves(g)]
    end
    return graphs_vac_vec
end

"""
    getGraphsG(graphs_vec::Vector{Vector{Graph}}) -> Vector{Vector{GraphG}}

Generate all unique graphs with external legs from regular graphs.

# Arguments
- `graphs_vec::Vector{Vector{Graph}}`: Nested vector of graphs by order

# Returns
- `Vector{Vector{GraphG}}`: Collection of GraphG objects by order

# Description
For each graph in the input collection, generates all possible ways to
attach external legs (marked vertices), filtering out duplicates and
ensuring the resulting graphs have no generalized leaves and no leaves
except at the marked vertices.

# Examples
```julia
# Generate graphs with external legs from regular graphs
leg_graphs = getGraphsG(all_graphs_by_order)
```

# Note
Results are cached to disk to avoid repeated computation.
"""
function getGraphsG(graphs_vec::Vector{Vector{Graph}})::Vector{Vector{GraphG}}
    
    graphsG_vec = [[GraphG(graphs_vec[1][1].g,[1,1])]]

    for n in eachindex(graphs_vec[2:end])
        g_vec = graphs_vec[2:end][n]
        gG_vec = GraphG[]
        fileName = data_dir()*"/GraphFiles/graphsG_"*string(n)*".jld2"

        ### load if available
        if isfile(fileName)
            println("load GraphsG for n="*string(n))
            gG_vec = load_object(fileName)
        
        ### or compute if not
        else
            println("compute GraphsG for n="*string(n))
            for g in g_vec
                ggG_vec = GraphG[]
                for jp in vertices(g.g), j in vertices(g.g)
                    if noLeavesExceptAt(g.g,[j,jp])
                        gG_test = GraphG(g.g,[j,jp])
                        if !hasGeneralizedLeaves(gG_test)
                            isNew = true
                            Threads.@threads for ggG in ggG_vec
                                if isIsomorph(ggG,gG_test) 
                                    isNew = false
                                    break
                                end
                            end 
                            if isNew push!(ggG_vec,gG_test) end
                        end
                    end
                end
                append!(gG_vec,ggG_vec)
            end
            save_object(fileName,gG_vec)
        end

        push!(graphsG_vec,gG_vec)
    end

    return graphsG_vec
end


"""
    merge_data_files()

Initialize and merge required data files if they don't exist.
"""
function merge_data_files()
    if !isfile(data_dir()*"/GraphFiles/graphs_12.jld2")
        println("Begin merging DynHTE datafiles:")
        println("This will only happen the very first time this function is used. It can take a while.")

        println("merging graphs12 ...")
        save_object(data_dir()*"/GraphFiles/graphs_12.jld2",
            vcat(load_object(data_dir()*"/GraphFiles/graphs_12a.jld2"),
                 load_object(data_dir()*"/GraphFiles/graphs_12b.jld2")))
    end

    if !isfile(data_dir()*"/GraphFiles/graphsG_12.jld2")
        println("merging graphsG12 ...")
        save_object(data_dir()*"/GraphFiles/graphsG_12.jld2",
            vcat(load_object(data_dir()*"/GraphFiles/graphsG_12a.jld2"),
                 load_object(data_dir()*"/GraphFiles/graphsG_12b.jld2"),
                 load_object(data_dir()*"/GraphFiles/graphsG_12c.jld2"),
                 load_object(data_dir()*"/GraphFiles/graphsG_12d.jld2")))
    end

    for sstring in ["S1half","S1"]
        if !isfile(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_11.jld2")
            println(sstring*": merging C_11 ...")
            save_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_11.jld2",
                vcat(load_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_11a.jld2"),
                     load_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_11b.jld2")))
        end

        if !isfile(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_12.jld2")
            println(sstring*": merging C_12 ...")
            save_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_12.jld2",
                vcat(load_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_12a.jld2"),
                     load_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_12b.jld2"),
                     load_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_12c.jld2"),
                     load_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_12d.jld2"),
                     load_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_12e.jld2"),
                     load_object(data_dir()*"/GraphEvaluations/Spin_"*sstring*"/C_12f.jld2")))
        end
    end
end


################################################################################
###### START ###################################################################
################################################################################

###### generation of basic graphs with two or fewer leaves
#graphs_vec = [ [Graph(SimpleWeightedGraph{Int64,Int64}(diagm([Int64(0)])))] ]      ## the single-vertex graph with 0 edges
#save_object(data_dir()*"/GraphFiles/graphs_0.jld2",graphs_vec[end])
#for nn in 1:8
#    @show nn
#    push!(graphs_vec,getAllGraphsNextOrder(graphs_vec[end]));       
#end
#gplot(graphs_vec[1+6])

##repeat the evaluation of this line to generate all graphs
#push!(graphs_vec,getAllGraphsNextOrder(graphs_vec[end]));  #graphs11 with 46384 graphs use 36.4MiB


###### if already generated: load graphs_1,2,3,...,nmax
#nmax=12 
#graphs_vec = [load_object(data_dir()*"/GraphFiles/graphs_"*string(n)*".jld2") for n in 0:nmax]

###### generate/load lists of vac-graphs and graphsG
#graphsVac_vec = getVacGraphs(graphs_vec)
#graphsG_vec = getGraphsG(graphs_vec)