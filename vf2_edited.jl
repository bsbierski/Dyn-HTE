#Copyright (c) 2015: Seth Bromberger and other contributors. Copyright (c) 2012: John Myles White and other contributors.

#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

using Graphs
import Graphs.Experimental.IsomorphismAlgorithm
import Graphs.Experimental.GraphMorphismProblem
import Graphs.Experimental.SubGraphIsomorphismProblem
import Graphs.Experimental.IsomorphismProblem
import Graphs.Experimental.could_have_isomorph
import Graphs.SimpleGraphs.AbstractSimpleGraph
"""
    VF2

An empty concrete type used to dispatch to [`vf2`](@ref) isomorphism functions.
"""
struct VF2 <: IsomorphismAlgorithm end




"""
    VF2State{G, T}

Structure that is internally used by vf2
"""
struct VF2State{G,T}
    g1::G
    g2::SimpleWeightedGraph{T,T}
    core_1::Vector{T}
    core_2::Vector{T}
    in_1::Vector{T}
    in_2::Vector{T}
    out_1::Vector{T}
    out_2::Vector{T}

    function VF2State(g1::G, g2::SimpleWeightedGraph{T,T}) where {G<:AbstractSimpleGraph{T}} where {T<:Integer}
        n1 = nv(g1)
        n2 = nv(g2)
        core_1 = zeros(T, n1)
        core_2 = zeros(T, n2)
        in_1 = zeros(T, n1)
        in_2 = zeros(T, n2)
        out_1 = zeros(T, n1)
        out_2 = zeros(T, n2)

        return new{G,T}(g1, g2, core_1, core_2, in_1, in_2, out_1, out_2)
    end
end

function identity_relation(i,j)
    return true
end


"""
    vf2(callback, g1, g2, problemtype; vertex_relation=nothing, edge_relation=nothing)

Iterate over all isomorphism between the graphs `g1` (or subgraphs thereof) and `g2`.
The problem that is solved depends on the value of `problemtype`:
- IsomorphismProblem(): Only isomorphisms between the whole graph `g1` and `g2` are considered.
- SubGraphIsomorphismProblem(): All isomorphism between subgraphs of `g1` and `g2` are considered.
- InducedSubGraphIsomorphismProblem(): All isomorphism between vertex induced subgraphs of `g1` and `g2` are considered.

Upon finding an isomorphism, the function `callback` is called with a vector `vmap` as an argument.
`vmap` is a vector where `vmap[v] == u` means that vertex `v` in `g2` is mapped to vertex `u` in `g1`.
If the algorithm should look for another isomorphism, then this function should return `true`.

### Optional Arguments
- `vertex_relation`: A binary function that takes a vertex from `g1` and one from `g2`. An
    isomorphism only exists if this function returns `true` for all matched vertices.
- `edge_relation`: A binary function that takes an edge from `g1` and one from `g2`. An
    isomorphism only exists if this function returns `true` for all matched edges.

### References
Luigi P. Cordella, Pasquale Foggia, Carlo Sansone, Mario Vento
“A (Sub)Graph Isomorphism Algorithm for Matching Large Graphs”
"""
function vf2(
    callback::Function,
    g1::G,
    g2::SimpleWeightedGraph,
    problemtype::GraphMorphismProblem;
    vertex_relation::F1 = identity_relation,
    edge_relation::F2 = identity_relation,
    jL1::Int8 = Int8(0),
    jL2::Int8 = Int8(0),
    jG1::Int8 = Int8(0),
    jG2::Int8 = Int8(0) 
    ) where {F1,F2,G<:AbstractSimpleGraph}
    if nv(g1) < nv(g2) || (problemtype == IsomorphismProblem() && nv(g1) != nv(g2))
        return nothing
    end

    state = VF2State(g1, g2)
    depth = 1

    if problemtype == SubGraphIsomorphismProblem()
        #the external vertices have to match either way.
        u = jL1
        v = jG1
        vf2update_state!(state, u, v, depth)
        depth = depth +1 

        if jL1 != jL2 
            u = jL2
            v = jG2
            vf2update_state!(state, u, v, depth)
            depth = depth +1 
        end

    end
 
    vf2match!(state, depth, callback, problemtype, vertex_relation, edge_relation)
    return nothing
end

"""
    vf2check_feasibility(u, v, state, problemtype, vertex_relation, edge_relation)

Check whether two vertices of G₁ and G₂ can be matched. Used by [`vf2match!`](@ref).
"""
function vf2check_feasibility(
    u,
    v,
    state::VF2State,
    problemtype,
    vertex_relation::F1,
    edge_relation::F2,
    ) where {F1,F2}
    @inline function vf2rule_pred(u, v, state::VF2State, problemtype)
        if problemtype != SubGraphIsomorphismProblem()
            @inbounds for u2 in inneighbors(state.g1, u)
                if state.core_1[u2] != 0
                    found = false
                    # TODO can probably be replaced with has_edge for better performance
                    for v2 in inneighbors(state.g2, v)
                        if state.core_1[u2] == v2
                            found = true
                            break
                        end
                    end
                    found || return false
                end
            end
        end
        @inbounds for v2 in inneighbors(state.g2, v)
            if state.core_2[v2] != 0
                found = false
                for u2 in inneighbors(state.g1, u)
                    if state.core_2[v2] == u2
                        found = true
                        break
                    end
                end
                found || return false
            end
        end
        return true
    end

    @inline function vf2rule_succ(u, v, state::VF2State, problemtype)
        if problemtype != SubGraphIsomorphismProblem()
            @inbounds for u2 in outneighbors(state.g1, u)
                if state.core_1[u2] != 0
                    found = false
                    for v2 in outneighbors(state.g2, v)
                        if state.core_1[u2] == v2
                            found = true
                            break
                        end
                    end
                    found || return false
                end
            end
        end
        found = false
        @inbounds for v2 in outneighbors(state.g2, v)
            if state.core_2[v2] != 0
                found = false
                for u2 in outneighbors(state.g1, u)
                    if state.core_2[v2] == u2
                        found = true
                        break
                    end
                end
                found || return false
            end
        end
        return true
    end

    @inline function vf2rule_in(u, v, state::VF2State, problemtype)
        count1 = 0
        count2 = 0
        @inbounds for u2 in outneighbors(state.g1, u)
            if state.in_1[u2] != 0 && state.core_1[u2] == 0
                count1 += 1
            end
        end
        @inbounds for v2 in outneighbors(state.g2, v)
            if state.in_2[v2] != 0 && state.core_2[v2] == 0
                count2 += 1
            end
        end
        if problemtype == IsomorphismProblem()
            count1 == count2 || return false
        else
            count1 >= count2 || return false
        end
        count1 = 0
        count2 = 0
        @inbounds for u2 in inneighbors(state.g1, u)
            if state.in_1[u2] != 0 && state.core_1[u2] == 0
                count1 += 1
            end
        end
        @inbounds for v2 in inneighbors(state.g2, v)
            if state.in_2[v2] != 0 && state.core_2[v2] == 0
                count2 += 1
            end
        end
        problemtype == IsomorphismProblem() && return count1 == count2

        return count1 >= count2
    end

    @inline function vf2rule_out(u, v, state::VF2State, problemtype)
        count1 = 0
        count2 = 0
        @inbounds for u2 in outneighbors(state.g1, u)
            if state.out_1[u2] != 0 && state.core_1[u2] == 0
                count1 += 1
            end
        end
        @inbounds for v2 in outneighbors(state.g2, v)
            if state.out_2[v2] != 0 && state.core_2[v2] == 0
                count2 += 1
            end
        end
        if problemtype == IsomorphismProblem()
            count1 == count2 || return false
        else
            count1 >= count2 || return false
        end

        count1 = 0
        count2 = 0
        @inbounds for u2 in inneighbors(state.g1, u)
            if state.out_1[u2] != 0 && state.core_1[u2] == 0
                count1 += 1
            end
        end
        @inbounds for v2 in inneighbors(state.g2, v)
            if state.out_2[v2] != 0 && state.core_2[v2] == 0
                count2 += 1
            end
        end
        problemtype == IsomorphismProblem() && return count1 == count2

        return count1 >= count2
    end

    @inline function vf2rule_new(u, v, state::VF2State, problemtype)
        problemtype == SubGraphIsomorphismProblem() && return true
        count1 = 0
        count2 = 0
        @inbounds for u2 in inneighbors(state.g1, u)
            if state.in_1[u2] == 0 && state.out_1[u2] == 0
                count1 += 1
            end
        end
        @inbounds for v2 in inneighbors(state.g2, v)
            if state.in_2[v2] == 0 && state.out_2[v2] == 0
                count2 += 1
            end
        end
        if problemtype == IsomorphismProblem()
            count1 == count2 || return false
        else
            count1 >= count2 || return false
        end
        count1 = 0
        count2 = 0
        @inbounds for u2 in outneighbors(state.g1, u)
            if state.in_1[u2] == 0 && state.out_1[u2] == 0
                count1 += 1
            end
        end
        @inbounds for v2 in outneighbors(state.g2, v)
            if state.in_2[v2] == 0 && state.out_2[v2] == 0
                count2 += 1
            end
        end
        problemtype == IsomorphismProblem() && return count1 == count2

        return count1 >= count2
    end

    @inline function vf2rule_self_loops(u, v, state, problemtype)
        u_selflooped = has_edge(state.g1, u, u)
        v_selflooped = has_edge(state.g2, v, v)

        if problemtype == SubGraphIsomorphismProblem()
            return u_selflooped || !v_selflooped
        end
        return u_selflooped == v_selflooped
    end

    syntactic_feasability =
        vf2rule_pred(u, v, state, problemtype) &&
        vf2rule_succ(u, v, state, problemtype) &&
        vf2rule_in(u, v, state, problemtype) &&
        vf2rule_out(u, v, state, problemtype) &&
        vf2rule_new(u, v, state, problemtype) &&
        vf2rule_self_loops(u, v, state, problemtype)
    syntactic_feasability || return false

    #if vertex_relation != nothing
        vertex_relation(u, v) || return false
    #end
    #if edge_relation != nothing
        E1 = edgetype(state.g1)
        E2 = edgetype(state.g2)
        for u2 in outneighbors(state.g1, u)
            state.core_1[u2] == 0 && continue
            v2 = state.core_1[u2]
            edge_relation(E1(u, u2), E2(v, v2)) || return false
        end
        for u2 in inneighbors(state.g1, u)
            state.core_1[u2] == 0 && continue
            v2 = state.core_1[u2]
            edge_relation(E1(u2, u), E2(v2, v)) || return false
        end
   # end
    return true
end

"""
    vf2update_state!(state, u, v, depth)

Update state before recursing. Helper function for [`vf2match!`](@ref).
"""
function vf2update_state!(state::VF2State, u::Int8, v::Int8, depth)::Nothing
    @inbounds begin
        state.core_1[u] = v
        state.core_2[v] = u
        for w in outneighbors(state.g1, u)
            if state.out_1[w] == 0
                state.out_1[w] = depth
            end
        end
        for w in inneighbors(state.g1, u)
            if state.in_1[w] == 0
                state.in_1[w] = depth
            end
        end
        for w in outneighbors(state.g2, v)
            if state.out_2[w] == 0
                state.out_2[w] = depth
            end
        end
        for w in inneighbors(state.g2, v)
            if state.in_2[w] == 0
                state.in_2[w] = depth
            end
        end
    end
end

"""
    vf2reset_state!(state, u, v, depth)

Reset state after returning from recursion. Helper function for [`vf2match!`](@ref).
"""
function vf2reset_state!(state::VF2State, u::Int8, v::Int8, depth)::Nothing
    @inbounds begin
        state.core_1[u] = 0
        state.core_2[v] = 0
        for w in outneighbors(state.g1, u)
            if state.out_1[w] == depth
                state.out_1[w] = 0
            end
        end
        for w in inneighbors(state.g1, u)
            if state.in_1[w] == depth
                state.in_1[w] = 0
            end
        end
        for w in outneighbors(state.g2, v)
            if state.out_2[w] == depth
                state.out_2[w] = 0
            end
        end
        for w in inneighbors(state.g2, v)
            if state.in_2[w] == depth
                state.in_2[w] = 0
            end
        end
    end
end

"""
    vf2match!(state, depth, callback, problemtype, vertex_relation, edge_relation)

Perform isomorphic subgraph matching. Called by [`vf2`](@ref).
"""
function vf2match!(
    state,
    depth,
    callback::Function,
    problemtype::GraphMorphismProblem,
    vertex_relation::F1,
    edge_relation::F2,
    ) where {F1,F2}
    n1 = Int8(nv(state.g1))
    n2 = Int8(nv(state.g2))
    # if all vertices of G₂ are matched we call the callback function. If the
    # algorithm should look for another isomorphism then callback has to return true
    if depth > n2
        keepgoing = callback(state.core_2)
        return keepgoing
    end
    # First we try if there is a pair of unmatched vertices u∈G₁ v∈G₂ that are connected
    # by an edge going out of the set M(s) of already matched vertices
    found_pair = false
    v = 0
    @inbounds for j in Int8.(1:n2)
        if state.out_2[j] != 0 && state.core_2[j] == 0
            v = j
            break
        end
    end
    if v != 0
        @inbounds for u in Int8.(1:n1)
            if state.out_1[u] != 0 && state.core_1[u] == 0
                found_pair = true
                if vf2check_feasibility(
                    u, v, state, problemtype, vertex_relation, edge_relation
                )
                    vf2update_state!(state, u, v, depth)
                    keepgoing = vf2match!(
                        state,
                        depth + 1,
                        callback,
                        problemtype,
                        vertex_relation,
                        edge_relation,
                    )
                    keepgoing || return false
                    vf2reset_state!(state, u, v, depth)
                end
            end
        end
    end
    found_pair && return true
    # If that is not the case we try if there is a pair of unmatched vertices u∈G₁ v∈G₂ that
    # are connected  by an edge coming in from the set M(s) of already matched vertices
    v = 0
    @inbounds for j in Int8.(1:n2)
        if state.in_2[j] != 0 && state.core_2[j] == 0
            v = j
            break
        end
    end
    if v != 0
        @inbounds for u in Int8.(1:n1)
            if state.in_1[u] != 0 && state.core_1[u] == 0
                found_pair = true
                if vf2check_feasibility(
                    u, v, state, problemtype, vertex_relation, edge_relation
                )
                    vf2update_state!(state, u, v, depth)
                    keepgoing = vf2match!(
                        state,
                        depth + 1,
                        callback,
                        problemtype,
                        vertex_relation,
                        edge_relation,
                    )
                    keepgoing || return false
                    vf2reset_state!(state, u, v, depth)
                end
            end
        end
    end
    found_pair && return true
    # If this is also not the case, we try all pairs of vertices u∈G₁ v∈G₂ that are not
    # yet matched
    v = 0
    @inbounds for j in Int8.(1:n2)
        if state.core_2[j] == 0
            v = j
            break
        end
    end
    if v != 0
        @inbounds for u in Int8.(1:n1)
            if state.core_1[u] == 0
                if vf2check_feasibility(
                    u, v, state, problemtype, vertex_relation, edge_relation
                )
                    vf2update_state!(state, u, v, depth)
                    keepgoing = vf2match!(
                        state,
                        depth + 1,
                        callback,
                        problemtype,
                        vertex_relation,
                        edge_relation,
                    )
                    keepgoing || return false
                    vf2reset_state!(state, u, v, depth)
                end
            end
        end
    end
    return true
end

function count_subgraphisomorph(
    g1::AbstractGraph,
    g2::AbstractGraph;
    vertex_relation::F1=identity_relation,
    edge_relation::F2=identity_relation,
    jL1::Int8,
    jL2::Int8,
    jG1::Int8,
    jG2::Int8
    )::Int where {F1,F2}
    

   #=  result = 0
    callback(vmap) = (result += 1; return true) =#

  
  result = Ref(0) 
   function callback2(result::Ref{Int})::Bool
      result[] += 1; 
      return true
   end

   function callback(vmap)::Bool
    return callback2(result)
 end
   
    vf2(
        callback,
        g1,
        g2,
        SubGraphIsomorphismProblem();
        vertex_relation=vertex_relation,
        edge_relation=edge_relation,
        jL1= jL1, jL2 = jL2,
        jG1 = jG1 ,jG2 = jG2
    )
    
    return result[]
    #return result
end


function has_isomorph(
    g1::AbstractGraph,
    g2::AbstractGraph;
    vertex_relation::F1=identity_relation,
    edge_relation::F2=identity_relation,
    )::Bool where {F1,F2}
    !could_have_isomorph(g1, g2) && return false

    result = false
    callback(vmap) = (result = true; return false)
    vf2(
        callback,
        g1,
        g2,
        IsomorphismProblem();
        vertex_relation=vertex_relation,
        edge_relation=edge_relation,
    )
    return result
end

