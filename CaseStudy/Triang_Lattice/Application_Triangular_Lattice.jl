using JLD2
#using RobustPad

include("../../Embedding.jl")
include("../../GraphGeneration.jl")
include("../../LatticeGraphs.jl")
include("../../ConvenienceFunctions.jl") 
#specify max order
max_order = 11

#LOAD FILES 
#-------------------------------------------------------------------------------------

#generate list of graphs
graphs_vec = [load_object("GraphFiles/graphs_"*string(nn)*".jld2") for nn in 0:max_order];
gG_vec = vcat(getGraphsG([graphs_vec[1]]), [load_object("GraphFiles/graphsG_"*string(nn)*".jld2") for nn in 1:max_order]);
## identify which gG have the same underlying simple-graph structure. Precalculate Symmetry factors. 
gG_vec_unique = give_unique_gG_vec(gG_vec);

#create vector of all lower order dictionaries
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,max_order+1) 
   
#load dictionaries of all lower orders C_Dict_vec 
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object("GraphEvaluations/Spin_S1half/C_"*string(ord)*".jld2")
end 
#-----------------------------------

#1. Define lattice ball for embedding (it is enough for embedding of max_order graphs to have ball radius L=max_order)
L = max_order
lattice,LatGraph,center_sites = getLattice_Ball(L,"triang");
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.1,fontsize=7,nodeshape=:rect,curves=false))

#2.Compute all correlations in the lattice
#@time Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);
@load "CaseStudy/Triang_Lattice/Correlation_Data_L11.jld2" Correlators

#check the susceptibility with 10.1103/PhysRevB.53.14228
(brillouin_zone_cut([(0.0,0.0) (0.0,0.0) ;(0.0,0.0)  (0.0,0.0)],Correlators,lattice,center_sites)[1]*4)[:,1].*[factorial(n)*4^n for n in 0:max_order]


#3. Fourier Transform

### Compute A 2D Brillouin zone cut: 
N = 100
kx = range(-0+0.001,2pi+0.001,length=N)
ky = range(-0+0.001,2pi+0.001,length=N) 
kmat = [[1,1/sqrt(3)].*x .+  [1,-1/sqrt(3)].*y for x in kx, y in ky ]
kmat = [(y,x) for x in kx, y in ky ]
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
x = -1.5
padetype = [5,6]
evaluate(y) = eval_correlator_LR_continuous_pad(y, x, padetype); #define evaluation function
struc = real.(evaluate.(structurefactor));
p = Plots.heatmap(kx,ky,struc,clims=(0,0.25))

### Compute A 2D Brillouin zone cut: 
N = 100
kx = range(-2pi,2pi,length=N)
ky = range(-2pi,2pi,length=N) #[0.] #for chains
kmat = [(1,1/sqrt(3)).*x .+  (1,-1/sqrt(3)).*y for x in kx, y in ky ]


### Calculate the full S(ω,k)
x = 1.5
N = 50
kx = range(-0+0.001,2pi+0.001,length=N)
ky = range(-0+0.001,2pi+0.001,length=N) 
w_vec = collect(-5:10/100:5.0)
kmat = [[1,1/sqrt(3)].*x .+  [1,-1/sqrt(3)].*y for x in kx, y in ky ]
Swk = zeros(N,N,length(w_vec));
Threads.@threads for i = 1:N
for  j = 1:N
    Swk[i,j,:] .=  get_JSkw_mat_finitex("total","pade",x,[kmat[i,j]],w_vec,0.02,1,2,200,false,Correlators,lattice,center_sites)[1,:]
end
end
4*sum(struc)*((maximum(w_vec)-minimum(w_vec))/length(w_vec))*(4*pi^2/N^2)/(4*pi^2)

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
w_ind = 51
Plots.heatmap(kx,ky,Swk[:,:,w_ind], title = L"$S(\omega,k)$ at $\omega$ ="*string(w_vec[w_ind]))


############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone
#---triangular
path = [(0+0.001,0+0.001),(0,2pi/(sqrt(3))),(2pi/(3),2pi/(sqrt(3))),(0+0.001,0+0.001)]
pathticks = ["Γ","M","K","Γ"]




#Generate the path 

Nk = 100
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);

###### S(k,w) heatmap
using CairoMakie

x = 1.5
w_vec = collect(0.01:0.0314/2:5.0)
JSkw_mat = get_JSkw_mat_finitex("total","pade",x,kvec,w_vec,0.02,1,2,200,false,Correlators,lattice,center_sites)


fig = Figure(fontsize=25,resolution = (900,500));
ax=Axis(fig[1,1],limits=(0,Nk+1,0,5),ylabel=L"\omega/J=w",title="Triangular Lattice: x="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25);
hm=CairoMakie.heatmap!(ax,[k for k in 1:Nk+1],w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,0.5),highclip=:white);
ax.xticks = (kticks_positioins,pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=40, label = L"J S(k,w)")
resize_to_layout!(fig);
display(fig)

#save("TriangX"*string(x)*".pdf",fig)

#Compare Static susceptibility to QMC
using CSV
using DataFrames
using Plots

# Load the dataset
df = CSV.File("CaseStudy/Triang_Lattice/TriangularDiagMC.csv") |> DataFrame


function extract_floats(col)
    return filter(!isnothing, [tryparse(Float64, x) for x in col if !ismissing(x)])
end

# Extract columns for plotting
x = [extract_floats(df[!, i]) for i in 1:2:7]
y = [extract_floats(df[!, i+1]) for i in 1:2:7]
labels = ["T = 0.5", "T = 1.0", "T = 0.375", "T = 2.0"]

# Create the plot
p = Plots.plot()
for i in 1:4
    Plots.plot!(p, Nk * x[i], y[i], seriestype = :scatter, label = labels[i], markersize = 2)
end
Plots.ylims!(p, (0, 1.9))
Plots.title!(p, "Triangular Static susceptibility")

####Plot the susceptibility along the high symmetry path for various temperatures 
#---triangular
path = [(0,0),(2pi/(3),2pi/(sqrt(3))),(0,2pi/(sqrt(3))),(0.,0.)]
pathticks = ["Γ","K","M","Γ"]

BrillPath[60][:,1]

#Generate the path 
Nk = 100;
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk);
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);

#Plot the path
Tvec = [2.0,1.0,0.5,0.375];
padetypes = [[5,6],[6,4],[4,6],[5,5]];
for padetype in padetypes
for T in Tvec
evaluate(x) = eval_correlator_LR_continuous_pad(x,-1/T,padetype);
cut = 1/T* real.(evaluate.(BrillPath));
Plots.plot!(p,0:Nk,cut,label=false, xticks=(kticks_positioins,pathticks)#= ,yrange = [0,2],xrange = [900,2800] =#)
end
end
p

Plots.plot!(p,[0,1],[-1,-1], label = "Pade", color=:black)
Plots.scatter!(p,[0,1],[-1,-1], label = "QMC", color=:black)


#Check Quantum to Classical Correspondence

# Define the tight-binding dispersion relation for the triangular lattice
function triangular_band(kx, ky, t=1.0, a=1.0)
    return -t * (cos(kx * a) + 2*cos(0.5 * kx * a) * cos(0.5 * sqrt(3) * ky * a))
end

data = [ -0.42 ./ (triangular_band(k[1], k[2]) .- 1.9) for k in kvec]
data = [ -0.44 ./ (triangular_band(k[1], k[2]) .- 1.79) for k in kvec]
q = Plots.plot(p,0:Nk,data, label = "Pade", color=:black)




