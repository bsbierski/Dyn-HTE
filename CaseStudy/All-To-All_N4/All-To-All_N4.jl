using JLD2,DelimitedFiles
include("../../plotConventions.jl")
include("../../Embedding.jl")
include("../../LatticeGraphs.jl")
include("../../ConvenienceFunctions.jl") 

### maximum expansion order and spin-length
n_max = 11
spin_string = "S1" #"S1half" #"S3half"

### load list of unique graphs
gG_vec_unique = give_unique_gG_vec(n_max);

### load dictionaries of all lower orders C_Dict_vec 
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,n_max+1);
for ord = 0:n_max
    C_Dict_vec[ord+1]  = load_object("GraphEvaluations/Spin_"*spin_string*"/C_"*string(ord)*".jld2")
end 

### define lattice for embedding
LatGraph = complete_graph(4)
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.11,fontsize=4,nodeshape=:rect,curves=false))

### compute all correlations in the lattice (or load them)
c11 = Calculate_Correlator_fast(LatGraph,1,1,n_max,gG_vec_unique,C_Dict_vec)
c12 = Calculate_Correlator_fast(LatGraph,1,2,n_max,gG_vec_unique,C_Dict_vec)
