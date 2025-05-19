### all tools revolving around abstract graphs and their visualization
using Random, Parameters, JLD2
using Graphs, SimpleWeightedGraphs
using SparseArrays, Combinatorics, LinearAlgebra
using GraphRecipes,Plots
plt_empty = plot(label="",axis=([], false))
graphsInRow = 6 #for plotting

include("Structs.jl")
include("vf2_edited.jl") 

####### various helper functions ################### 
function load_dyn_hte_graphs(spin_length::Number,max_order::Int)::Dyn_HTE_Graphs
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
        C_Dict_vec[ord+1]  = load_object("GraphEvaluations/"*S_string*"/C_"*string(ord)*".jld2")
    end 
    Dyn_HTE_Graphs(S,gG_vec_unique,C_Dict_vec)
end


function split_vec(vec::Vector,part::Int,parts::Int)
    """ splits vector in parts (keep longtail), returns the chunk and its start and end indices """
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

function degeneracy(g::SimpleWeightedGraph{Int64, Int64})::Int
    """ get degeneracy factor prod_b(m_b!) for bonds b in SimpleWeightedGraph """
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

function isIsomorph(g1::Graph,g2::Graph)::Bool
    """ for graphs: check if g1 ~ g2 """
    ### convert g1,2 to SimpleGraphs
    g1_simple = toSimpleGraph(g1.g)
    g2_simple = toSimpleGraph(g2.g)

    ### prepare isomomorphism check on g1,2_simple respecting edge weigths in g1,g2 
    edge_relation(b1,b2) = (g1.g.weights[src(b1),dst(b1)] == g2.g.weights[src(b2),dst(b2)])
    
    return Graphs.Experimental.has_isomorph(g1_simple,g2_simple,edge_relation=edge_relation)
end
function isIsomorph(gG1::GraphG,gG2::GraphG)::Bool
    """ for graphsG: check if gG1 ~ gG2 with (j1,j1')=(j2,j2') or (j1,j1')=(j2',j2)"""
    
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

function is_symmetric(gG::GraphG)::Bool
    """
    find if gG is a symmetric graph with respect to switching the two external legs (fix the mapping of the two external vertices)
    """
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

function findg(g::Graph,g_vec::Vector{Graph})::Int
    """ search g in g_vec, if found return index, else return 0"""
    for k in eachindex(g_vec)
        if isIsomorph(g,g_vec[k])
            return k
        end
    end
    return 0
end
function findg(gG::GraphG,gG_vec::Vector{GraphG})::Int
    """ find gG in gG_vec, if found return index, else return 0"""
    for k in eachindex(gG_vec)
        if isIsomorph(gG,gG_vec[k])
            return k
        end
    end
    return 0
end

function addOneEdge(g::Graph)::Vector{Graph}
    """ take Graph g and add one edge in all possible ways (strenghten existing edge or connect to new vertex) """
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

function totalEdges(g::SimpleWeightedGraph{Int64, Int64})::Int
    return Int(sum(g.weights)/2)
end

function removeVerticesWithoutEdge!(g::Graph)
    """ do not affect trivial graph with single vertex """
    for j in 2:nv(g)
        if length(outneighbors(g,j))==0
            rem_vertex!(g,j)
        end
    end
end

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

function splitToConnectedComp(g::SimpleWeightedGraph{Int64, Int64})::Vector{Graph}
    """ split vacuum graph into connected components, drop isolated vertices (this is used after subtraction) """
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
function hasGeneralizedLeaves(g::Graph)::Bool
    """ check if a Graph g has generalized leaves (used for vacuum graphs)"""
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
function hasGeneralizedLeaves(gG::GraphG)::Bool
    """ check if a GraphG gG has generalized leaves """
    
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

function symmetryFactor(g::Graph)::Int
    """ symmetry Factor of graph, this is the number auf automorphisms respecting edge weights"""
    gg = copy(g.g)

    ### convert gg to SimpleGraphs
    gg_simple = toSimpleGraph(gg)

    ### prepare auto-momorphism check on gg_simple respecting edge weigths in gg 
    edge_relation(b1,b2) = (gg.weights[src(b1),dst(b1)] == gg.weights[src(b2),dst(b2)])

    return Graphs.Experimental.count_isomorph(gg_simple,gg_simple,edge_relation=edge_relation)
end
function symmetryFactor(gG::GraphG)::Int
    """ symmetry Factor of graphG, this is the number of graph automorphisms respecting edge weights and not touching external indices"""
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

####### graph plotting ##################################################################
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

###### graph generation ############################################################
function getAllGraphsNextOrder(g_vec::Vector{Graph})::Vector{Graph}
    """ iterate to find all connected multi-graphs of one order higher than those provided, skip graphs with more than two leaves """
    
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
    save_object("GraphFiles/graphs_"*string(totalEdges(g_new_vec[end].g))*".jld2",g_new_vec)
    return g_new_vec
end

function getVacGraphs(graphs_vec::Vector{Vector{Graph}})::Vector{Vector{Graph}}
    graphs_vac_vec = copy(graphs_vec)
    for m in 1:length(graphs_vac_vec)-1
        graphs_vac_vec[m+1] = [g for g in graphs_vac_vec[m+1] if !hasGeneralizedLeaves(g)]
    end
    return graphs_vac_vec
end

function getGraphsG(graphs_vec::Vector{Vector{Graph}})::Vector{Vector{GraphG}}
    """ load graphsG or compute GraphsG from graphs_vec, need to remove equivalent ways of adding terminals """
    graphsG_vec = [[GraphG(graphs_vec[1][1].g,[1,1])]]

    for n in eachindex(graphs_vec[2:end])
        g_vec = graphs_vec[2:end][n]
        gG_vec = GraphG[]
        fileName = "GraphFiles/graphsG_"*string(n)*".jld2"

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

### if GraphFiles/graphs_12.jld2 has not yet been merged from its <100Mb parts a,b, then merge and save it
if !isfile("GraphFiles/graphs_12.jld2")
    println("merging graphs12 ...")
    save_object("GraphFiles/graphs_12.jld2",vcat(load_object("GraphFiles/graphs_12a.jld2"),load_object("GraphFiles/graphs_12b.jld2")))
end
### if GraphFiles/graphsG_12.jld2 has not yet been merged from its <100Mb parts a,b,c,d then merge and save it
if !isfile("GraphFiles/graphsG_12.jld2")
    println("merging graphsG12 ...")
    save_object("GraphFiles/graphsG_12.jld2",vcat(load_object("GraphFiles/graphsG_12a.jld2"),load_object("GraphFiles/graphsG_12b.jld2"),load_object("GraphFiles/graphsG_12c.jld2"),load_object("GraphFiles/graphsG_12d.jld2")))
end

### if GraphEvaluations C_11.jld2 and C_12.jld2 do not yet exist, merge it from its parts
for sstring in ["S1half","S1"]

    if !isfile("GraphEvaluations/Spin_"*sstring*"/C_11.jld2")
        println(sstring*": merging C_11 ...")
        save_object("GraphEvaluations/Spin_"*sstring*"/C_11.jld2",vcat(     load_object("GraphEvaluations/Spin_"*sstring*"/C_11a.jld2"),
                                                                            load_object("GraphEvaluations/Spin_"*sstring*"/C_11b.jld2")            
                                                                    ))
    end

    if !isfile("GraphEvaluations/Spin_"*sstring*"/C_12.jld2")
        println(sstring*": merging C_12 ...")
        save_object("GraphEvaluations/Spin_"*sstring*"/C_12.jld2",vcat( load_object("GraphEvaluations/Spin_"*sstring*"/C_12a.jld2"),
                                                                        load_object("GraphEvaluations/Spin_"*sstring*"/C_12b.jld2"),
                                                                        load_object("GraphEvaluations/Spin_"*sstring*"/C_12c.jld2"),
                                                                        load_object("GraphEvaluations/Spin_"*sstring*"/C_12d.jld2"),
                                                                        load_object("GraphEvaluations/Spin_"*sstring*"/C_12e.jld2"),
                                                                        load_object("GraphEvaluations/Spin_"*sstring*"/C_12f.jld2")            
                                                                    ))
    end

end


################################################################################
###### START ###################################################################
################################################################################

###### generation of basic graphs with two or fewer leaves
#graphs_vec = [ [Graph(SimpleWeightedGraph{Int64,Int64}(diagm([Int64(0)])))] ]      ## the single-vertex graph with 0 edges
#save_object("GraphFiles/graphs_0.jld2",graphs_vec[end])
#for nn in 1:8
#    @show nn
#    push!(graphs_vec,getAllGraphsNextOrder(graphs_vec[end]));       
#end
#gplot(graphs_vec[1+6])

##repeat the evaluation of this line to generate all graphs
#push!(graphs_vec,getAllGraphsNextOrder(graphs_vec[end]));  #graphs11 with 46384 graphs use 36.4MiB


###### if already generated: load graphs_1,2,3,...,nmax
nmax=12 
graphs_vec = [load_object("GraphFiles/graphs_"*string(n)*".jld2") for n in 0:nmax]

###### generate/load lists of vac-graphs and graphsG
graphsVac_vec = getVacGraphs(graphs_vec)
graphsG_vec = getGraphsG(graphs_vec)




function whichParts(parts::Int,vec,miss::Vector{Int})::Vector{Int}
    res = Int[]
    for part in 1:parts
        @show part
        chunk,a,b = split_vec(vec,part,parts)
        for m in miss
            if m in collect(a:b)
                append!(res,part)
                break
            end
        end
    end
    return res
end