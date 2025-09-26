### all tools revolving around abstract graphs and their visualization
plt_empty = plot(label="", axis=([], false))
graphsInRow = 6 #for plotting

####### various helper functions ################### 
"""
    load_dyn_hte_graphs(spin_length::Number, max_order::Int; verbose = false) -> Dyn_HTE_Graphs

Load pre-computed graph-evaluation data for given spin-length (1/2 or 1) and up to order max_order. Work with unique graphs.

# Arguments
- `spin_length::Number`: The spin length (e.g., 1/2 or 1)
- `max_order::Int`: Maximum expansion order to load
- `verbose::Bool=false`: Whether to print loading information

# Returns
- `Dyn_HTE_Graphs`: Structure containing unique graphs and their coefficients
"""
function load_dyn_hte_graphs(spin_length::Number, max_order::Int; verbose=false)::Dyn_HTE_Graphs
    #first check if all required data files are available, if not merge them
    merge_data_files()

    S = rationalize(spin_length)
    if S == 1 // 2
        S_string = "Spin_S1half"
    elseif S == 1 // 1
        S_string = "Spin_S1"
    else
        throw(error("Spin-length " * string(spin_length) * " is not yet implemented"))
    end

    #load list of unique graphs
    gG_vec_unique = give_unique_gG_vec(max_order)

    #create vector of all lower order dictionaries
    C_Dict_vec = Vector{Vector{Vector{Rational{Int128}}}}(undef, max_order + 1)
    #load dictionaries of all lower orders C_Dict_vec 
    for ord = 0:max_order
        if verbose
            println("loading C_" * string(ord) * " for S=" * S_string)
        end
        C_Dict_vec[ord+1] = load_object(data_dir() * "/GraphEvaluations/" * S_string * "/C_" * string(ord) * ".jld2")
    end
    Dyn_HTE_Graphs(S, gG_vec_unique, C_Dict_vec)
end


"""
    degeneracy(g::SimpleWeightedGraph{Int64, Int64}) -> Int

Calculate the degeneracy factor of a weighted graph (product of the factorials of edge multiplicities)
"""
function degeneracy(g::SimpleWeightedGraph{Int64,Int64})::Int
    mat = g.weights
    deg = 1
    n = nv(g)
    for row in 1:n
        for col in row+1:n
            deg *= factorial(mat[row, col])
        end
    end
    return deg
end

"""
    isIsomorph(g1::Graph, g2::Graph) -> Bool

Check if two Graph objects are isomorphic, respecting edge weights.
"""
function isIsomorph(g1::Graph, g2::Graph)::Bool

    ### convert g1,2 to SimpleGraphs
    g1_simple = toSimpleGraph(g1.g)
    g2_simple = toSimpleGraph(g2.g)

    ### prepare isomomorphism check on g1,2_simple respecting edge weigths in g1,g2 
    edge_relation(b1, b2) = (g1.g.weights[src(b1), dst(b1)] == g2.g.weights[src(b2), dst(b2)])

    return Graphs.Experimental.has_isomorph(g1_simple, g2_simple, edge_relation=edge_relation)
end

"""
    isIsomorph(gG1::GraphG, gG2::GraphG) -> Bool

Check if two GraphG objects (graphs with terminal vertices) are isomorphic. Terminal vertices can match j1=j1', j2=j2' or flipped j1=j2', j2=j1'
"""
function isIsomorph(gG1::GraphG, gG2::GraphG)::Bool


    ### add terminal vertices via edge with weight 100 
    gG1_ext = copy(gG1.g)
    gG2_ext = copy(gG2.g)

    for j in 1:2
        add_vertex!(gG1_ext)
        add_edge!(gG1_ext, gG1.jjp[j], nv(gG1_ext), 100)

        add_vertex!(gG2_ext)
        add_edge!(gG2_ext, gG2.jjp[j], nv(gG2_ext), 100)
    end

    ### convert g1G,gG2 to SimpleGraphs
    gG1_ext_simple = toSimpleGraph(gG1_ext)
    gG2_ext_simple = toSimpleGraph(gG2_ext)

    ### prepare isomomorphism check on g1,2_simple respecting edge weigths in g1,g2 
    edge_relation(b1, b2) = (gG1_ext.weights[src(b1), dst(b1)] == gG2_ext.weights[src(b2), dst(b2)])

    return Graphs.Experimental.has_isomorph(gG1_ext_simple, gG2_ext_simple, edge_relation=edge_relation)
end

"""
    is_symmetric(gG::GraphG) -> Bool

Check if gG is a symmetric graph with respect to switching the two external vertices.
"""
function is_symmetric(gG::GraphG)::Bool

    gg = gG.g
    gg_simple = toSimpleGraph(gg)

    edge_relation(b1, b2) = (gg.weights[src(b1), dst(b1)] == gg.weights[src(b2), dst(b2)])

    # finds if there is an isomorphism by only permuting the internal vertices between the graph and the graph with its external vertices flipped.
    count = count_subgraphisomorph(gg_simple, gg, edge_relation=edge_relation, jL1=gG.jjp[2], jL2=gG.jjp[1], jG1=gG.jjp[1], jG2=gG.jjp[2])

    if count > 0
        return true
    else
        return false
    end
end

"""
    totalEdges(g::SimpleWeightedGraph{Int64, Int64}) -> Int

Count the total number of edges (including multiplities) in a weighted graph.
"""
function totalEdges(g::SimpleWeightedGraph{Int64,Int64})::Int
    return Int(sum(g.weights) / 2)
end

"""
    toSimpleGraph(g::SimpleWeightedGraph{Int64, Int64}) -> SimpleGraph

Convert a weighted graph to a simple (unweighted) graph.
"""
function toSimpleGraph(g::SimpleWeightedGraph{Int64,Int64})::SimpleGraph
    ### drop the weights of a SimpleWeightedGraph to get a SimpleGraph
    g_simple = SimpleGraph(nv(g))
    for v1 in 1:nv(g), v2 in v1+1:nv(g)
        if g.weights[v1, v2] > 0
            add_edge!(g_simple, v1, v2)
        end
    end
    return g_simple
end


"""
    noLeavesExceptAt(g::SimpleWeightedGraph{Int64, Int64}, j_vec::Vector{Int64}=Int64[]) -> Bool

Check if a graph has leaves only at specified vertices in j_vec (at these vertices, there are external operators)
"""
function noLeavesExceptAt(g::SimpleWeightedGraph{Int64,Int64}, j_vec::Vector{Int64}=Int64[])::Bool
    ### check if a Graph g has leaves, exclude sites in j_vec from checking
    for i in vertices(g)
        bonds_i = sum(g.weights[i, :])
        if bonds_i < 2 && !(i in j_vec)
            return false
        end
    end
    return true
end

"""
    hasGeneralizedLeaves(g::Graph) -> Bool

Check if a graph has generalized leaves. A generalized leaf is a bridge edge (weight 1) whose removal would disconnect
the graph. This function checks for the existence of such edges.
"""
function hasGeneralizedLeaves(g::Graph)::Bool

    g_ext = copy(g.g)

    ### move through all edges e, if weigth=1 remove it and check if there is still just one connected component 
    for e in edges(g_ext)
        if get_weight(g_ext, src(e), dst(e)) == 1
            g_ext_rem = copy(g_ext)
            rem_edge!(g_ext_rem, e)
            g_ext_rem = toSimpleGraph(g_ext_rem)
            con_comp = connected_components(g_ext_rem)
            if length(con_comp) > 1
                return true
            end
        end
    end
    return false
end

"""
    hasGeneralizedLeaves(gG::GraphG) -> Bool

Check if a graphG with external legs has generalized leaves. A generalized leaf is a bridge edge (weight 1)
whose removal would create a separate component that still connects the
external vertices. This function checks for the existence of such edges.
"""
function hasGeneralizedLeaves(gG::GraphG)::Bool
    ### add terminal vertices via edge with weight 100 
    gG_ext = copy(gG.g)

    for j in 1:2
        add_vertex!(gG_ext)
        add_edge!(gG_ext, gG.jjp[j], nv(gG_ext), 100)
    end

    ### move through all edges e, if weigth=1 remove it, if j,jp are still connected see if we get a second disconnected component 
    for e in edges(gG_ext)
        if get_weight(gG_ext, src(e), dst(e)) == 1
            gG_ext_rem = copy(gG_ext)
            rem_edge!(gG_ext_rem, e)
            gG_ext_rem = toSimpleGraph(gG_ext_rem)
            con_comp = connected_components(gG_ext_rem)
            if length(con_comp) > 1
                Nv = nv(gG.g)
                for vs in con_comp
                    if issubset([Nv + 1, Nv + 2], vs)
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
The symmetry factor is the number of graph isomorphisms that respect the edge-multiplicities.
"""
function symmetryFactor(g::Graph)::Int

    gg = copy(g.g)

    ### convert gg to SimpleGraphs
    gg_simple = toSimpleGraph(gg)

    ### prepare auto-momorphism check on gg_simple respecting edge weigths in gg 
    edge_relation(b1, b2) = (gg.weights[src(b1), dst(b1)] == gg.weights[src(b2), dst(b2)])

    return Graphs.Experimental.count_isomorph(gg_simple, gg_simple, edge_relation=edge_relation)
end

"""
    symmetryFactor(gG::GraphG) -> Int

Calculate the symmetry factor of a graph with external legs.
The symmetry factor is the number of graph isomorphisms that keep the terminal vertices invariant and respect the edge-multiplicities.
"""
function symmetryFactor(gG::GraphG)::Int

    gg = copy(gG.g)

    ### add to gg two vertices at terminals j,j' with bond-weight 100,101
    add_vertex!(gg)
    add_edge!(gg, gG.jjp[1], nv(gg), 100)
    add_vertex!(gg)
    add_edge!(gg, gG.jjp[2], nv(gg), 101)

    ### convert gg,gg_flip to SimpleGraphs
    gg_simple = toSimpleGraph(gg)

    ### prepare auto-momorphism check on gg_simple respecting edge weigths in gg 
    edge_relation(b1, b2) = (gg.weights[src(b1), dst(b1)] == gg.weights[src(b2), dst(b2)])

    return Graphs.Experimental.count_isomorph(gg_simple, gg_simple, edge_relation=edge_relation)
end

"""
    gplot(gG::GraphG; title::String="", save::Bool=false) -> Plot
    gplot(g::Graph; title::String="", save::Bool=false) -> Plot

Create a visualization of a Graph or GraphG (graph with external legs, highlighted with red flags). Also takes vectors.
"""
function gplot(g::Graph; title::String="", save::Bool=false)
    Random.seed!(2011)
    m = collect(g.g.weights)
    n = size(m)[1]

    ### catch case of no bonds --- does not work yet!!!
    if maximum(m) == 0
        fig = graphplot(g, names=1:n, markersize=0.22, fontsize=12, nodeshape=:circle)
    else
        ### need to convert in multigraph language that yields the destination sites
        gg = Vector{Int}[]
        for src in 1:n-1
            dests = Int[]
            for dest in src+1:n
                if m[src, dest] > 0
                    append!(dests, [dest for _ in 1:m[src, dest]])
                end
            end
            append!(gg, [dests])
        end
        fig = graphplot(gg, names=[string(n) for n in 1:n], markercolor=:lightblue, markersize=0.15, fontsize=14, nodeshape=:circle, arrow=false, edgewidth=(s, d, w) -> 3, title=title)
    end

    ### save fig as png if filename supplied
    if save
        savefig("figs/" * title * ".png")
    end
    return fig
end
function gplot(gG::GraphG; title::String="", save::Bool=false)
    Random.seed!(2011)
    gplot = copy(gG.g) #the SimpleWeightedGraph for plotting with external legs

    ### add external bonds in red
    for j in gG.jjp
        add_vertex!(gplot)
        add_edge!(gplot, j, nv(gplot))
    end

    m = collect(gplot.weights)
    n = size(m)[1]

    function edgewidth_function(s, d, w)
        if d >= n - 1
            return 1
        end
        return 3
    end

    ### need to convert in multigraph language that yields the destination sites
    gg = Vector{Int}[]
    for src in 1:n-1
        dests = Int[]
        for dest in src+1:n
            if m[src, dest] > 0
                append!(dests, [dest for _ in 1:m[src, dest]])
            end
        end
        append!(gg, [dests])
    end
    fig = graphplot(gg, names=append!([string(n) for n in 1:n-2], ["j", "j'"]), markercolor=append!([:lightblue for _ in 1:n-2], [:red, :red]), nodesize=0.15, nodeweights=append!(ones(n - 2), [0.4, 0.4]), fontsize=14, nodeshape=append!([:circle for _ in 1:n-2], [:rect, :rect]), ew=edgewidth_function, arrow=false, title=title)

    ### save fig as png if filename supplied
    if save
        savefig("figs/" * title * ".png")
    end
    return fig
end
function gplot(g_vec::Vector{}; subtitle_vec::Vector{String}=["#" * string(pos) for pos in eachindex(g_vec)], title::String="", save::Bool=false)
    fig_vec = []
    for g_pos in 1:length(g_vec)
        fig = gplot(g_vec[g_pos]; title=subtitle_vec[g_pos])
        append!(fig_vec, [fig])
    end
    rows = Int(floor(length(fig_vec) / graphsInRow))
    if rows * graphsInRow < length(g_vec)
        rows += 1
    end
    fig = plot(fig_vec..., layout=(rows, graphsInRow), dpi=300, size=(300 * graphsInRow, 300 * rows), plot_title=title)
    if save
        savefig("figs/" * title * ".png")
    end
    return fig
end

"""
    getVacGraphs(graphs_vec::Vector{Vector{Graph}}) -> Vector{Vector{Graph}}

Filter a collection of graphs to remove those with generalized leaves.
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

Generate all graphs with external legs from regular graphs..
"""
function getGraphsG(graphs_vec::Vector{Vector{Graph}})::Vector{Vector{GraphG}}

    graphsG_vec = [[GraphG(graphs_vec[1][1].g, [1, 1])]]

    for n in eachindex(graphs_vec[2:end])
        g_vec = graphs_vec[2:end][n]
        gG_vec = GraphG[]
        fileName = data_dir() * "/GraphFiles/graphsG_" * string(n) * ".jld2"

        ### load if available
        if isfile(fileName)
            println("load GraphsG for n=" * string(n))
            gG_vec = load_object(fileName)

            ### or compute if not
        else
            println("compute GraphsG for n=" * string(n))
            for g in g_vec
                ggG_vec = GraphG[]
                for jp in vertices(g.g), j in vertices(g.g)
                    if noLeavesExceptAt(g.g, [j, jp])
                        gG_test = GraphG(g.g, [j, jp])
                        if !hasGeneralizedLeaves(gG_test)
                            isNew = true
                            Threads.@threads for ggG in ggG_vec
                                if isIsomorph(ggG, gG_test)
                                    isNew = false
                                    break
                                end
                            end
                            if isNew
                                push!(ggG_vec, gG_test)
                            end
                        end
                    end
                end
                append!(gG_vec, ggG_vec)
            end
            save_object(fileName, gG_vec)
        end

        push!(graphsG_vec, gG_vec)
    end

    return graphsG_vec
end

"""
    merge_data_files()

Initialize and merge required data files if they don't exist.
"""
function merge_data_files()
    if !isfile(data_dir() * "/GraphFiles/graphs_12.jld2")
        println("Begin merging DynHTE datafiles:")
        println("This will only happen the very first time this function is used. It can take a while.")

        println("merging graphs12 ...")
        save_object(data_dir() * "/GraphFiles/graphs_12.jld2",
            vcat(load_object(data_dir() * "/GraphFiles/graphs_12a.jld2"),
                load_object(data_dir() * "/GraphFiles/graphs_12b.jld2")))
    end

    if !isfile(data_dir() * "/GraphFiles/graphsG_12.jld2")
        println("merging graphsG12 ...")
        save_object(data_dir() * "/GraphFiles/graphsG_12.jld2",
            vcat(load_object(data_dir() * "/GraphFiles/graphsG_12a.jld2"),
                load_object(data_dir() * "/GraphFiles/graphsG_12b.jld2"),
                load_object(data_dir() * "/GraphFiles/graphsG_12c.jld2"),
                load_object(data_dir() * "/GraphFiles/graphsG_12d.jld2")))
    end

    for sstring in ["S1half", "S1"]
        if !isfile(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_11.jld2")
            println(sstring * ": merging C_11 ...")
            save_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_11.jld2",
                vcat(load_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_11a.jld2"),
                    load_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_11b.jld2")))
        end

        if !isfile(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_12.jld2")
            println(sstring * ": merging C_12 ...")
            save_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_12.jld2",
                vcat(load_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_12a.jld2"),
                    load_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_12b.jld2"),
                    load_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_12c.jld2"),
                    load_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_12d.jld2"),
                    load_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_12e.jld2"),
                    load_object(data_dir() * "/GraphEvaluations/Spin_" * sstring * "/C_12f.jld2")))
        end
    end
end
