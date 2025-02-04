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

#generate list of graphs
graphs_vec = [load_object("GraphFiles/graphs_"*string(nn)*".jld2") for nn in 0:max_order];
gG_vec = getGraphsG(graphs_vec);
## identify which gG have the same underlying simple-graph structure. Precalculate Symmetry factors. 
gG_vec_unique = give_unique_gG_vec(gG_vec);

#create vector of all lower order dictionaries
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,max_order+1) 
   
#load dictionaries of all lower orders C_Dict_vec 
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object("GraphFiles/GraphG_Lists/C_"*string(ord)*".jld2")
end 
#-----------------------------------

#1. Define lattice ball for embedding (it is enough for embedding of max_order graphs to have ball radius L=max_order)
L = max_order
lattice,LatGraph,center_sites = getLattice_Ball(L,"square");
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.1,fontsize=7,nodeshape=:rect,curves=false))

#2.Compute all correlations in the lattice
Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);



#3. Fourier Transform

### Compute A 2D Brillouin zone cut: 
N = 100
kx = range(-4pi,4pi,length=N)
ky = range(-4pi,4pi,length=N) #[0.] #for chains
kmat = [(y,x) for x in kx, y in ky ]
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
x = 1.2
padetype = [4,4]
evaluate(y) = eval_correlator_LR_continuous_pad(y, x, padetype); #define evaluation function
struc = real.(evaluate.(structurefactor));
p = Plots.heatmap(kx,ky,struc,clims=(0,1.))

############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone
#---square
path = [(0,0),(pi,0),(pi,pi),(0,0)]
pathticks = ["Γ","X","M","Γ"]

#---square path Sherman
path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
pathticks = ["K","X","M","K","Γ","X"]


#Generate the path 
Nk = 100
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);

###### S(k,w) heatmap
using CairoMakie

x = 2.0
w_vec = collect(0.01:0.0314/2:5.0)
JSkw_mat = get_JSkw_mat_finitex("total","pade",x,kvec,w_vec,0.02,0,3,200,false,Correlators,lattice,center_sites)


fig = Figure(fontsize=25,resolution = (900,500));
ax=Axis(fig[1,1],limits=(0,Nk+1,0,5),ylabel=L"\omega/J=w",title="Square Lattice: x="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25);
hm=CairoMakie.heatmap!(ax,[k for k in 1:Nk+1],w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,1.0),highclip=:white);
ax.xticks = (kticks_positioins,pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=40, label = L"J S(k,w)")
resize_to_layout!(fig);
display(fig)

#Add T=0 QMC Data and Linear Spin Wave Theory
@load "CaseStudy/Square_Lattice/QMC_square_T0_data.jld2" QMC_axis QMC_data
disp = [sqrt(1-1/4*(cos(k[1])+cos(k[2]))^2)  for k in kvec]

CairoMakie.plot!(ax,QMC_axis*(Nk+0.5),1*QMC_data, color = :pink, alpha = 0.45,label = L"T=0 \text{QMC-Gap}")
CairoMakie.plot!(ax,[k for k in 1:Nk+1],2.4*disp, color = :orange, alpha = 0.45,label = L"T=0 \text{LSWT-Dispersion}")

axislegend(ax)
display(fig)

#save("SquareX"*string(x)*".pdf",fig)

