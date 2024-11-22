### all tools revolving around abstract graphs and their visualization
using Random, Parameters, JLD2
using Graphs, SimpleWeightedGraphs
using SparseArrays, Combinatorics, LinearAlgebra
using GraphRecipes,Plots
plt_empty = plot(label="",axis=([], false))
graphsInRow = 6 #for plotting

include("Structs.jl")

####### various helper functions ################### 
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
    """ for graphsG: check if gG1 ~ gG2 (j1,j1')=(j2,j2') or =(j2',j2)"""
    
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

function isSymmetric(gG::GraphG)::Bool
    """ check if gG is symmetric under exchange of j,j' """
    if gG.jjp[1]==gG.jjp[2] return true
    else
        gg = copy(gG.g)
        gg_flip = copy(gG.g)
        
        ### add to gg terminals j,j' and j',j with bond-weight 100,200
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
        
        ### prepare isomorphism check on gg <--> ggflip respecting the edge-weights
        edge_relation(b1,b2) = (gg.weights[src(b1),dst(b1)] == gg_flip.weights[src(b2),dst(b2)])
        
        ### check if gG,j,j' is isomorph to gG,j',j
        return Graphs.Experimental.has_isomorph(gg_simple,gg_flip_simple,edge_relation=edge_relation)
    end
end

function findg(g::Graph,g_vec::Vector{Graph})::Int
    """ find g in g_vec, if found return index, else return 0"""
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
function noLeavesExceptAt(g::SimpleWeightedGraph{Int64, Int64},j_vec::Vector{Int}=Int[])::Bool
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

    ### add to gg two vertices at terminals j,j' with bond-weight 100,200
    add_vertex!(gg)
    add_edge!(gg,gG.jjp[1],nv(gg),100)
    add_vertex!(gg)
    add_edge!(gg,gG.jjp[2],nv(gg),200)

    ### convert gg,gg_flip to SimpleGraphs
    gg_simple = toSimpleGraph(gg)

    ### prepare auto-momorphism check on gg_simple respecting edge weigths in gg 
    edge_relation(b1,b2) = (gg.weights[src(b1),dst(b1)] == gg.weights[src(b2),dst(b2)])

    return Graphs.Experimental.count_isomorph(gg_simple,gg_simple,edge_relation=edge_relation)
end

####### graph plotting ##################################################################
function gplot(g::Graph;title::String="",save::Bool=false)
    Random.seed!(2017);
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
    Random.seed!(2017);
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
    for m in 1:length(graphs_vac_vec)
        graphs_vac_vec[m] = [g for g in graphs_vac_vec[m] if !hasGeneralizedLeaves(g)]
    end
    return graphs_vac_vec
end

function getGraphsG(graphs_vec::Vector{Vector{Graph}})::Vector{Vector{GraphG}}
    """ get GraphsG from graphs_vec, need to remove equivalent ways of adding terminals """
    graphsG_vec = [[GraphG(graphs_vec[1][1].g,[1,1])]]

    for n in eachindex(graphs_vec[2:end])
        g_vec = graphs_vec[2:end][n]
        gG_vec = GraphG[]

        println("get GraphsG for n="*string(n))
        
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
        
        push!(graphsG_vec,gG_vec)
    end

    return graphsG_vec
end

###### V-con subtraction (edge based!) ########################################
function getVconSubtractions(gG::GraphG,graphsG_vec::Vector{Vector{GraphG}},graphsVac_vec::Vector{Vector{Graph}})::Vector{VconSub}
    """ get V-con subtraction by keeping subsets of edges of gG to form gG_sub """
    ### obtain vector of all edges (for multi-edges this contains multiples)
    edge_vec = Vector{Int}[]
    for j = 1:nv(gG.g)
        for i in j+1:nv(gG.g)
            for _ in 1:gG.g.weights[i,j]
                push!(edge_vec,[i,j])
            end
        end
    end

    ### prepare results vec
    res_vec = []
    
    ### iterate true subset V' of bonds, subtract them from gc and put them to g
    for Vp in Iterators.drop(unique(powerset(edge_vec)),1)
        # initialize sub-graphs: gG as parent, g with all edges removes
        gG_SWG_sub,g_SWG_sub = copy(gG.g),copy(gG.g)
        for e in edges(g_SWG_sub)
            rem_edge!(g_SWG_sub,e)
        end

        ### for each edge b from V': move b from gc to g ### TODO: speed-up here by checking that j-j' still remains connected
        for b in Vp
            if !add_edge!(gG_SWG_sub,b...,get_weight(gG_SWG_sub,b...)-1) 
                rem_edge!(gG_SWG_sub,b...)
            end
            add_edge!(g_SWG_sub ,b...,get_weight(g_SWG_sub ,b...)+1)
        end

        ### gG_SWG_sub: keep only connected components which are not isolated vertices or isolated and connected to gG.jjp 
        conComp_gG = [v for v in connected_components(gG_SWG_sub) if (length(v)>1 || issubset(gG.jjp,v))]

        ### require: only one connected component also connected to terminals is left 
        gG_sub = nothing
        if length(conComp_gG)==1 && issubset(gG.jjp,conComp_gG[1])
            vs = conComp_gG[1]
            gG_SWG_sub_part = SimpleWeightedGraph(gG_SWG_sub.weights[vs,vs])
            gG_sub = GraphG(gG_SWG_sub_part,[findfirst(x->x==gG.jjp[1],vs),findfirst(x->x==gG.jjp[2],vs)])
        end
        
        ### proceed only if gG_sub is connected & has no generalized leaves and is thus isomorphic to one in list
        if !isnothing(gG_sub) && !hasGeneralizedLeaves(gG_sub)
            p = findg(gG_sub,graphsG_vec[1+totalEdges(gG_sub.g)])
            @assert(p>0)

            ###### treatment of product of vacuum graphs (it it is)
            ### split to connected component and prepare last entry of VconSub  
            g_sub_vec_candidate = splitToConnectedComp(g_SWG_sub)
            gVac_sub_vec = Tuple{Int,Int}[]
                    
            ### continue only if all g in g_sub_vec_candidate have no generalized leaves
            if maximum([hasGeneralizedLeaves(g) for g in g_sub_vec_candidate])==0
                ### find the isomorphic vacuum graphs
                for gVac_sub in g_sub_vec_candidate
                    m = totalEdges(gVac_sub.g)
                    for k in 1:length(graphsVac_vec[m+1])
                        if isIsomorph(gVac_sub,graphsVac_vec[m+1][k])
                            push!(gVac_sub_vec,(m,k))
                            break
                        end
                        if k==length(graphsVac_vec[m+1]) 
                            display(gplot(gVac_sub))
                            error("vacuum graph not found!") 
                        end
                    end
                end
            end

            ### push subtraction to res_vec
            if length(gVac_sub_vec)>0 
                res = [ (totalEdges(gG_sub.g),p) , gVac_sub_vec]
                push!(res_vec,res)
            end

        end
    end
    
    ### find multiplicity factors and return
    return [VconSub(count(x->x==sub,res_vec),sub[1],sub[2]) for sub in unique(res_vec)]
end

function plotVconSubtractions(gG::GraphG,subs::Vector{VconSub},graphsG_vec::Vector{Vector{GraphG}},graphsVac_vec::Vector{Vector{Graph}};save::Bool=false)
    
    if length(subs)==0
        maxNumberOfVacuumGraphs = 0
    else
        maxNumberOfVacuumGraphs = maximum([length(sub.g_sub_vec) for sub in subs])
    end
    plt_vec = append!([gplot(gG,title="subtract from:")],[gplot(graphsG_vec[sub.gG_sub[1]+1][sub.gG_sub[2]],title="factor="*string(sub.factor)) for sub in subs])

    for g_sub_pos in 1:maxNumberOfVacuumGraphs
        append!(plt_vec,[plt_empty])
        for sub in subs
            if length(sub.g_sub_vec)>=g_sub_pos
                append!(plt_vec,[ gplot(graphsVac_vec[sub.g_sub_vec[g_sub_pos][1]+1][sub.g_sub_vec[g_sub_pos][2]]) ] )
            else
                append!(plt_vec,[plt_empty])
            end
        end
    end

    fig = plot(plt_vec..., layout=(1+maxNumberOfVacuumGraphs,1+length(subs))
    ,dpi=300,size=(400*(1+length(subs)),400*(1+maxNumberOfVacuumGraphs)),plot_title="V-con subtraction")
    if save
        savefig("figs/"*"VConSub.png")
    end
    return fig
end

################################################################################
###### START ###################################################################
################################################################################

###### generation of basic graphs with two or fewer leaves
#nmax=7
#graphs_vec = [ [Graph(SimpleWeightedGraph(diagm([0])))] ]      ## the single-vertex graph with 0 edges
#save_object("GraphFiles/graphs_0.jld2",graphs_vec[end])
#for _ in 1:nmax
#    push!(graphs_vec,getAllGraphsNextOrder(graphs_vec[end]));       ##repeat the evaluation of this line to generate all graphs
#end
#gplot(graphs_vec[1+nmax])

###### take apart graphs_12 into a,b to comply with 100MB limit of git
#graphs_12 = load_object("GraphFiles/graphs_12.jld2")
#save_object("GraphFiles/graphs_12a.jld2",graphs_12[1:100000])
#save_object("GraphFiles/graphs_12b.jld2",graphs_12[100001:end])

###### if already generated: load graphs_1,2,3,...,10 or up to 12
#graphs_vec = [load_object("GraphFiles/graphs_"*string(n)*".jld2") for n in 0:10]
#graphs_vec = push!([load_object("GraphFiles/graphs_"*string(n)*".jld2") for n in 0:11],append!(load_object("GraphFiles/graphs_12a.jld2"),load_object("GraphFiles/graphs_12b.jld2")))

###### generate lists of vac-graphs and graphsG
#graphsVac_vec = getVacGraphs(graphs_vec)
#graphsG_vec = getGraphsG(graphs_vec)

###### plot low order graphs and some dedicated figs for paper
if false
    for nn in 2:5
        gplot(graphsG_vec[1+nn],title="graphsG_"*string(nn),save=true)
        gplot(graphsVac_vec[1+nn],title="graphsVac_"*string(nn),save=true)
    end
    gplot(union(graphsG_vec[0+1],graphsG_vec[1+1],graphsG_vec[2+1],graphsG_vec[3+1]),save=true,title="",subtitle_vec=["(n=0,r=1)","(1,1)","(2,1)","(2,2)","(2,3)","(3,1)","(3,2)","(3,3)","(3,4)","(3,5)","(3,6)","(3,7)"])
end

####### test the V-con subtraction with graphG [example from manuscript with gG labeled by (6,90)]
if false
    #gG = graphsG_vec[4+1][6]
    #gG = graphsG_vec[7+1][156]
    gG = graphsG_vec[6+1][90]
    #gG = graphsG_vec[5+1][36]

    gplot(gG)
    subs=getVconSubtractions(gG,graphsG_vec,graphsVac_vec)
    plotVconSubtractions(gG,subs,graphsG_vec,graphsVac_vec,save=true)
end


###### test tetramer g12 embeddings and subtractions
if false
    i,ip=1,1
    L = getLattice(4,"all-to-all")
    display(graphplot(L,names=1:nv(L),markersize=0.14,fontsize=10,nodeshape=:rect,curves=false))

    es=[e(L,i,ip,graphsG_vec[1+5][r]) for r in 1:length(graphsG_vec[1+5])]'
    ds=[degeneracy(graphsG_vec[1+5][r].g) for r in 1:length(graphsG_vec[1+5])]'

    gTest_vec = [graphsG_vec[1+5][r] for r in 1:length(graphsG_vec[1+5]) if e(L,i,ip,graphsG_vec[1+5][r])>0]
    display(gplot(gTest_vec;subtitle_vec=(["e="*string(es[r])*", d="*string(ds[r]) for r in 1:length(graphsG_vec[1+5]) if es[r]>0]),title="N=4 cycle"))

    #C_5=load_object("GraphFiles/GraphG_Lists/C_5.jld2")
    #[C_5[r] for r in 1:length(gG_vec[1+5]) if e(L,1,2,gG_vec[1+5][r])>0]

    #for gG in gTest_vec
    #    subs=getVconSubtractions(gG,graphsG_vec,graphsVac_vec)
    #    display(plotVconSubtractions(gG,subs,graphsG_vec,graphsVac_vec))
    #end
end