using JLD2,DelimitedFiles,Graphs,Plots
using DynHTE

### maximum expansion order and spin-length
n_max = 6
spin_length = 1

### define lattice for embedding (here: all-to-all with N=4 sites, any other lattice graph can be built)
Lgraph = complete_graph(4)
display(graphplot(Lgraph,names=1:nv(Lgraph),markersize=0.11,fontsize=4,nodeshape=:rect,curves=false))

### load graphG evaluations 
hte_graphs = load_dyn_hte_graphs(spin_length,n_max);

### compute all correlations in the lattice 
c_iipDyn_mat = get_c_iipDyn_mat(Lgraph,[1],hte_graphs)

### test if uniform susceptibility is purely static (m=0 only)
c_iipDyn_mat[1,1]+c_iipDyn_mat[2,1]+c_iipDyn_mat[3,1]+c_iipDyn_mat[4,1]


### test Dyn-HTE of local and non-local correlator (c.f. Eq. (A3) and (A4) in PRB-paper )
c_iipDyn_mat[1,1]
c_iipDyn_mat[2,1]