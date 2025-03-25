using JLD2
#using RobustPad

include("../../Embedding.jl")
include("../../GraphGeneration.jl")
include("../../LatticeGraphs.jl")
include("../../ConvenienceFunctions.jl") 
#specify max order
max_order = 8

#LOAD FILES 
#-------------------------------------------------------------------------------------
#load list of unique graphs
gG_vec_unique = give_unique_gG_vec(max_order);

#create vector of all lower order dictionaries
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,max_order+1) ;
#load dictionaries of all lower orders C_Dict_vec 
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object("GraphEvaluations/Spin_S1half/C_"*string(ord)*".jld2")
end 
#-----------------------------------

#1. Define lattice ball for embedding (it is enough for embedding of max_order graphs to have ball radius L=max_order)
L = max_order
lattice,LatGraph,center_sites = getLattice_Ball(L,"pyrochlore");
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))



#2.Compute all correlations in the lattice
@time Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);


#check the susceptibility with 10.1103/PhysRevB.53.14228
(brillouin_zone_cut([(0.0,0.0,0.0) (0.0,0.0,0.0) ;(0.0,0.0,0.0)  (0.0,0.0,0.0)],Correlators,lattice,center_sites)[1]*4)[:,1].*[factorial(n)*4^n for n in 0:max_order]



### Compute A 2D Brillouin zone cut: 
N = 100
kx = range(-8pi,8pi,length=N)
ky = range(-8pi,8pi,length=N) #[0.] #for chains
kmat = [(x,x,y) for x in kx, y in ky ]
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
x = -5
padetype = [4,4]
evaluate(y) = eval_correlator_LR_continuous_pad_exp(y, x, padetype,3); #define evaluation function
struc = ( evaluate.(structurefactor));
p = Plots.heatmap(kx,ky,abs(x)*transpose(struc); clims=(0,4),aspect_ratio=1/sqrt(3) , xlims = [-8pi,8pi])




############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone
#---pyrochlore (2,2, l)
path = [(4pi,4pi,-4pi+0.01),(4pi,4pi,0+0.01),(4pi,4pi,4pi+0.01)]
pathticks = ["Γ","W","Γ"]

#---pyrochlore (h,h, 2)
path = [(-4pi,-4pi,4pi+0.01),(0,0,4pi+0.01),(4pi,4pi,4pi+0.01)]
pathticks = ["Γ","W","Γ"]




#Generate the path 

Nk = 100
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);



###### S(k,w) heatmap
using CairoMakie

x = 3.
w_vec = collect(-5:0.0614/2:5.0)
JSkw_mat = get_JSkw_mat_finitex("total","pade",x,kvec,w_vec,0.02,1,2,200,false,Correlators,lattice,center_sites)


fig = Figure(fontsize=25,resolution = (900,500));
ax=Axis(fig[1,1],limits=(0,Nk+1,0,5),ylabel=L"\omega/J=w",title="Pyrochlore): x="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25);
hm=CairoMakie.heatmap!(ax,[k for k in 1:Nk+1],w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,1),highclip=:white);
ax.xticks = (kticks_positioins,pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=40, label = L"J S(k,w)")
resize_to_layout!(fig);
display(fig)

