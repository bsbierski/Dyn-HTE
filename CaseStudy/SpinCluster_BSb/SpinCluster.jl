using JLD2,DelimitedFiles
include("../../plotConventions.jl")
include("../../LatticeGraphs.jl")
include("../../Embedding.jl")
include("../../ConvenienceFunctions.jl") 

### maximum expansion order and spin-length
n_max = 12
spin_length = 1

### define lattice for embedding (here: all-to-all with N=4 sites)
Lgraph = complete_graph(4)
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.11,fontsize=4,nodeshape=:rect,curves=false))

### load graphG evaluations 
hte_graphs = load_dyn_hte_graphs(spin_length,n_max);

if false ### plot the C_n vs graphG index
    plt=plot(yscale=:log10,xlabel="graph index",ylabel="denom C_n",title="spin length S=$spin_length")
    for n in n_max:-1:8
        tmp = hte_graphs.c_dict[n+1]
        #tmp = load_object("GraphEvaluations/Spin_"*spin_string*"/Vac_"*string(n)*".jld2")
        plot!(plt,[abs(denominator(c[1]))+1 for c in tmp],label="n=$n")
    end
    display(plt)
    #savefig(plt,"C_n_S"*string(spin_length)*".png")
end

### compute all correlations in the lattice 
c_iipDyn_mat = get_c_iipDyn_mat(Lgraph,[1],hte_graphs)

### test if uniform susceptibility is pureyl static (m=0 only)
c_iipDyn_mat[1,1]+c_iipDyn_mat[2,1]+c_iipDyn_mat[3,1]+c_iipDyn_mat[4,1]+c_iipDyn_mat[5,1]+c_iipDyn_mat[6,1]+c_iipDyn_mat[7,1]+c_iipDyn_mat[8,1]

c_iipDyn_mat[1,1]
c_iipDyn_mat[2,1]
