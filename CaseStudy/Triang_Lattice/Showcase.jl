


# order of perturbation theory
L = 4
spin_length = 1/2

# import graphs
hte_graphs = load_dyn_hte_graphs(spin_length,L);

# import lattice
hte_lattice = getLattice(L,"square");

# plot lattice
display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

@time c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs);


N = 100
kx = range(-2pi,2pi,length=N)
ky = range(-2pi,2pi,length=N) 
kmat = [(y,x) for x in kx, y in ky ];
c_kDyn_mat =  get_c_kDyn(kmat,c_iipDyn_mat,hte_lattice);

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
x = -1.0
padetype = [2,2]
evaluate(y) = eval_correlator_LR_continuous_pad(y, x, padetype); #define evaluation function
struc = real.(evaluate.(c_kDyn_mat));
p = Plots.heatmap(kx,ky,struc,clims=(0,0.4))