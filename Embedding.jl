### find embedding factors
include("Structs.jl")
include("GraphGeneration.jl")
include("vf2_edited.jl") 

function is_symmetric(gG::GraphG)::Bool
    """
    find if gG is a symmetric graph with respect to switching the two external legs
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

function is_simple_isomorphic(gG1::GraphG,gG2::GraphG)::Bool
    """
    check if the underlying simple graphs of gG1 and gG2 are isomorphic
    """
 
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
#@save "GraphFiles/unique_gG_vec_0.jld2" vector
function give_unique_gG_vec(gG_vec::Vector{Vector{GraphG}})
    """
    gives unique_gG_vec with the structure
        [maxorder,
        [ 
        [gG,[gG_index_1,gG_index_2,...], dist ],
        ...
        ]
        where the gG_index_1,2,3 identify the GraphG of the same simple graph structure as gG. These indices have the structure
        gG_index_1 = [order + 1 ,index,symmetryfactor,is_symmetric] 
        dist = graph distance between external legs of gG
    """
    
    maxorder = length(gG_vec) - 1

    # try to load the file. if it does not exist try to load the file of one less order
    file_path = "GraphFiles/unique_gG_vec_$maxorder.jld2"

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
    @save "GraphFiles/unique_gG_vec_$maxorder.jld2" unique_gG_vec

    return unique_gG_vec

end


function give_unique_gG_vec(max_order::Int)
    """
    loads unique_gG_vec
    """

    # try to load the file. if it does not exist try to load the file of one less order
    file_path = "GraphFiles/unique_gG_vec_$max_order"*".jld2"

    if isfile(file_path)
        unique_gG_vec = load_object(file_path) 
        return unique_gG_vec
    else
        if max_order == 12
                ##Combine the 4 parts for order 12
                graphlist = Vector{Vector{unique_Graph}}(undef,4)
                for part = 1:4
                    @load "GraphFiles/unique_gG_vec_$max_order"*"_$part"*".jld2" unique_graphs_12_part
                    graphlist[part] = unique_Graphs(12,unique_graphs_12_part).graphs
                end
                    combined = vcat(graphlist...) 
                    combined_unique = unique_Graphs(max_order,combined)
                    @save "GraphFiles/unique_gG_vec_12.jld2" combined_unique
            return combined_unique
        end
        throw(ArgumentError("No unique graph file available for order $max_order"))
    end

end


function e_fast(LL::SimpleGraph{Int},j::Int,jp::Int,gG::GraphG)::Int
    """
    find embedding factor e 
    - lattice L (needs to be chosen large enough to avoid boundary effects!)
    - lattice site indices jjp=[j,j'] can be [i,i'] (or [i',i] if gG is not symmetric under exchange of i <--> i')
    - embedding of GraphG gG

    assumes that the distance j-jp is smaller or equal to the distance of external vertices of gG
    """

    numSubIsos = count_subgraphisomorph(LL,gG.g,jL1 = j,jL2 = jp,jG1 = gG.jjp[1],jG2 = gG.jjp[2])
    return numSubIsos 
end

function Calculate_Correlator_fast(L::SimpleGraph{Int},ext_j1::Int,ext_j2::Int,max_order::Int,gG_vec_unique::unique_Graphs,C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}})::Vector{Vector{Rational{Int128}}}
    """ Calculate the coefficients of (-x)^n for TG_ii'(iν_m) from embedding factors of only the unique simple graphs and the gG's symmetry factors """    

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

       # println("$index th graph embeding factor = $emb_fac")
        

        #### now we sum overall graphG eqivalent to the unique Gg
        for graph in unique_Gg.gG_vec
            g_order = graph.order #order
            gG_vec_index = graph.index #index
            symmetry_factor = graph.symmetry_factor#symmetry factor
            is_symmetric = graph.is_symmetric  #bool if the graph is symmetric
           
            fac = 2
            if  is_symmetric
                fac = 1
            end

            #look up the value of the graph from C_Dict_vec
            look_up_dict = C_Dict_vec[g_order+1][gG_vec_index]

            
            
            result_array[g_order+1] .+= look_up_dict.*Int128(emb_fac/symmetry_factor*fac)
        end
    end

    return result_array
end


###### LEGACY FUNCTIONS
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

function Calculate_Correlator(L::SimpleGraph{Int},ext_j1::Int,ext_j2::Int,max_order,gG_vec::Vector{Vector{GraphG}},C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}})::Vector{Vector{Rational{Int128}}}
    """OLD FUNCTION: use Calculate_Correlator_fast"""

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
