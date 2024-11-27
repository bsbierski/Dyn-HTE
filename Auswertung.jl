using JLD2
using RobustPade

#import Pkg 
#Pkg.activate(@__DIR__) #activates the environment in the folder of the current file

include("Embedding.jl")
include("GraphGeneration.jl")
include("LatticeGraphs.jl")
include("ConvenienceFunctions.jl") 
#specify max order
max_order = 10

#LOAD FILES 
#-------------------------------------------------------------------------------------

#generate list of graphs
graphs_vec = [load_object("GraphFiles/graphs_"*string(nn)*".jld2") for nn in 0:max_order];
gG_vec = getGraphsG(graphs_vec);
## identify which gG have the same underlying simple Graphs structure. Precalculate Symmetry factors. 
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
@time Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);


#3. Fourier Transform

### Compute A 2D Brillouin zone cut: 
N = 100
kx = range(-2pi,2pi,length=N)
ky = range(-2pi,2pi,length=N) #[0.] #for chains
kmat = [(x,y) for x in kx, y in ky ]
structurefactor = brillouin_zone_cut(kmat,Correlators,lattice,center_sites);


### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
iomega = 0
JoverT = 1.5
padetype = [4,4]
evaluate(x) = eval_correlator_LR_continuous_pad_Mats(x,iomega, JoverT, padetype); #define evaluation function
struc = real.(evaluate.(structurefactor));
 p = heatmap(kx,ky,struc)

############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone

#---chain
path = [(0,0),(2pi,0)]
pathticks = ["Γ","Γ"]

#---triangular
path = [(0,0),(0,2pi/(sqrt(3))),(2pi/(3),1.0*2pi/(sqrt(3))),(0.,0.)]
pathticks = ["Γ","M","K","Γ"]

#---square
path = [(0,0),(pi,0),(pi,pi),(0,0)]
pathticks = ["Γ","X","M","Γ"]


#Generate the path 
Nk = 200
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);



####Plot the susceptibility along the high symmetry path for various temperatures 
Tvec = [2.0,1.0,0.5,0.375]
padetype = [4,4]

p= plot();
for T in Tvec
evaluate(x) = eval_correlator_LR_continuous_pad_Mats(x,0,1/T,padetype);
cut = 1/T* real.(evaluate.(BrillPath));
plot!(p,1:length(kvec),cut, label = "T = $T", xticks=(kticks_positioins,pathticks)#= ,yrange = [0,2],xrange = [900,2800] =#)
end
p


####Plot the Spectrum along the high symmetry path for a set temperature  
ωvec = -5:0.1:5
JoverT = 0.5
padetype = [4,4]
spectrum = transpose([imag(eval_correlator_LR_continuous_pad(x,ω,JoverT,padetype)) for x in BrillPath, ω in ωvec]);
heatmap(1:length(kvec),ωvec,spectrum, xticks=(kticks_positioins,pathticks),colormap = :RdBu )

