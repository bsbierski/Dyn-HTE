using Graphs, SimpleWeightedGraphs,Parameters

###### graph structures, all Graphs are: connected (and thus free of vertices without any edges)
struct Graph ### connected vacuum graph
    g::SimpleWeightedGraph{Int64, Int64}
    function Graph(g)
        @assert is_connected(g)
        new(g)
    end
end

mutable struct GraphG ### V-connected graphs for correlator G with TWO external legs at sites j,jp
    g::SimpleWeightedGraph{Int64, Int64}
    jjp::Vector{Int64}
    function GraphG(g,jjp)

        @assert is_connected(g)
        @assert nv(g) >= maximum(jjp)
        @assert length(jjp)==2

        new(g,jjp)
    end
end

@with_kw mutable struct VconSub ### single subtraction to evaluate V-connected GraphG
    factor::Int = 1
    gG_sub::Tuple{Int,Int}
    g_sub_vec::Vector{Tuple{Int,Int}}
end



struct GraphGinfo
    gG::GraphG                      # graphG
    jjp::Tuple{Int8,Int8}           # jjp
    num_v::Int8                     # number of vertices
    n::Int8                         # order of graph = number of all edges
    m::Int8                         # different edges
    edge_vec::Vector{Vector{Int8}}  # all different edges in their order as in gG
    mult_vec::Vector{Int8}          # multiplicity of edges
    v_size::Vector{Int8}            # lengths of operator string for each vertex
    praw_vec::Vector{Int8}          # raw ordering vec of edge operators [1,2,2,3]
    
    function GraphGinfo(gG::GraphG)
        """ from gG::GraphG to GraphGinfo """
        num_v = nv(gG.g)        
        n = totalEdges(gG.g)    
        m = length(edges(gG.g))

        ### operator ordering info
        edge_vec = Vector{Int}[]    
        mult_vec = Int[]            
        praw_vec = Int[]            # m operators ordered as the numbers in p_vec: 1=Ve1, 2=Ve2,...
        
        multi_edge_vec = collect(edges(gG.g))
        for pos in eachindex(multi_edge_vec)
            e = multi_edge_vec[pos]
            push!(edge_vec, [src(e),dst(e)])
            push!(mult_vec, get_weight(gG.g,src(e),dst(e)))
            for _ in 1:get_weight(gG.g,src(e),dst(e))
                push!(praw_vec,pos)
            end
        end
        
        ### find lengths of all vertex operator strings
        v_size=zeros(Int,num_v)
        for e_pos in eachindex(edge_vec)
            e = 1*edge_vec[e_pos]
            v_size[e[1]] += mult_vec[e_pos]
            v_size[e[2]] += mult_vec[e_pos]
        end
        v_size[gG.jjp[1]] +=1
        v_size[gG.jjp[2]] +=1

        ### build struct
        new(gG,(gG.jjp[1],gG.jjp[2]),num_v,n,m,edge_vec,mult_vec,v_size,praw_vec)
    end
end

struct Struct_γ
    """ structure for all flavors to sum over """
    meγ::Vector{Vector{Int8}}   # multi-edge-operator resolved flavors, e.g. [[1,-1],[0],[1,0,-1]]
    vγ::Vector{Vector{Int8}}    # vertex resolved flavors (same as in Config_pγv, but only in initial order)
    Πx_fac::Int8                # multiplicity under Πx
end

mutable struct Config_pγv
    """ struct to keep track of operator permutations, edge flavors, γ-vecs for vertices, and their trace-evals """
    vγ::Vector{Vector{Int8}}    # γ-vecs for vertices
    pv::Vector{Vector{Int8}}    # vertex resolved operator order [3,5,8,9] with 0 the external vertex 
    vtr::Vector{Int}            # traces for vγ
    pos_ext::Vector{Int8}       # position of the external operators in vertex string order
    prod::Int                   # full trace
end

