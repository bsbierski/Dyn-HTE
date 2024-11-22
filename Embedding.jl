# find embedding factors
using Random, Parameters
using Graphs, SimpleWeightedGraphs,  SparseArrays, Combinatorics, LinearAlgebra
using GraphRecipes,Plots

include("Structs.jl")
include("GraphGeneration.jl")
include("vf2_edited.jl") 



function e(L::SimpleGraph{Int},j::Int,jp::Int,gG::GraphG)::Int
    """
    OLD FUNCTION use e_fast 
    find embedding factor e 
    - lattice L (needs to be chosen large enough to avoid boundary effects)
    - lattice site indices jjp=[j,j'] can be [i,i'] (or [i',i] if gG is not symmetric under exchange of i <--> i')
    - embedding of GraphG gG
    """
 
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




function is_symmetric(gG::GraphG)::Bool
    """
    find if gG1 is a symmetric graph with respect to switching the two external legs
    """
    gg1 = gG.g
    gg_simple = toSimpleGraph(gg1)
   

    edge_relation(b1,b2) = (gg1.weights[src(b1),dst(b1)] == gg1.weights[src(b2),dst(b2)])

    # finds if there is an isomorphism by only permuting the internal vertices between the graph and the graph with its external vertices flipped.
    count = count_subgraphisomorph(gg_simple,gg1,edge_relation=edge_relation,jL1 = gG.jjp[2],jL2 = gG.jjp[1],jG1 = gG.jjp[1],jG2 = gG.jjp[2])
    
    if count >0 
        return true
    else
        return false

    end
end


function is_isomorphic(gG1::GraphG,gG2::GraphG)::Bool
    """
    find if the underlying simple graphs of gG1 and gG2 are isomorphic
    """
 
    ### convert gg,gg_flip to SimpleGraphs
    gg1_simple = toSimpleGraph(gG1.g)
    

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


        return has_isomorph(gg1_simple,gG2.g ; vertex_relation= vertex_relation)
end

###Initialize the File if it does not exist yet.
#vector = [0,[[gG_vec[1][1],[[0,1,1,true]],0]]]
#@save "GraphFiles/unique_gG_vec_0.jld2" vector
function give_unique_gG_vec(gG_vec)
    """
    gives unique_gG_vec with the structurefactor
        [maxorder,
        [ 
        [gG,[gG_index_1,gG_index_2,...], dist ],
        ...
        ]
        where the indices have the structurefactor
        gG_index_1 = [order + 1 ,index,symmetryfactor,is_symmetric] 
        dist = distance between external legs
    """
    
    maxorder = length(gG_vec) - 1

    # try to load the file. if it does not exist try to load the file of one less order
    file_path = "GraphFiles/unique_gG_vec_$maxorder.jld2"

    if isfile(file_path)
        unique_gG_vec = load_object(file_path) 
        return unique_gG_vec
    else
        # if the order has not been calculated: try one order less
        unique_gG_vec =   give_unique_gG_vec(gG_vec[1:end-1])
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

        
        # find first isomorphic graph to gg that is already in the unique list. There is at most one. 
        unique_index = findfirst(x->is_isomorphic(gg,x[1]), unique_gG_vec[2])
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
    @save "GraphFiles/unique_gG_vec_$maxorder.jld2" unique_gG_vec

    return unique_gG_vec

end

function e_fast(LL::SimpleGraph{Int64},j::Int,jp::Int,gG::GraphG,issymmetric::Bool)::Int
    """
    find embedding factor e 
    - lattice L (needs to be chosen large enough to avoid boundary effects!)
    - lattice site indices jjp=[j,j'] can be [i,i'] (or [i',i] if gG is not symmetric under exchange of i <--> i')
    - embedding of GraphG gG

    ##assumes that the distance j-jp is smaller or equal to the distance of external vertices of gG
    """
    fac = 2
  if  issymmetric
    fac = 1
   end

    numSubIsos = count_subgraphisomorph(LL,gG.g,jL1 = j,jL2 = jp,jG1 = gG.jjp[1],jG2 = gG.jjp[2])
    return numSubIsos*fac 
end