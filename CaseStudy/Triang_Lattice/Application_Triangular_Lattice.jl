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
lattice,LatGraph,center_sites = getLattice_Ball(L,"triang");
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

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


#Generate the path 
p = plot()
Nk = 100;
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk);
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);
BrillPath[38][:,1]
#Plot the path
Tvec = [2.0,1.0,0.5,0.375];
padetypes = [[5,5]];
for padetype in padetypes
for T in Tvec
evaluate(x) = eval_correlator_LR_continuous_pad(x,-1/T,padetype);
cut = 1/T* real.(evaluate.(BrillPath));
Plots.plot!(p,0:Nk,cut,label=false, xticks=(kticks_positioins,pathticks)#= ,yrange = [0,2],xrange = [900,2800] =#)
end
end
p

p = plot()
Tvec = [2.0,1.0,0.5,0.375,0.1];
padetypes = [[5,5]];
for padetype in padetypes
for T in Tvec
evaluate(x) = eval_correlator_LR_continuous_pad_tanh(x,-1/T,padetype,4);
cut =  1/T* real.(evaluate.(BrillPath));
Plots.plot!(p,1:Nk+1,cut,label=false, xticks=(kticks_positioins,pathticks) , linestyle = :dash #= ,yrange = [0,2],xrange = [900,2800] =#)
end
end
p


for T in Tvec
    vals = [1/T*GTriang(1/T, k[1], k[2]) for k = kvec]
    Plots.plot!(p,1:Nk+1,vals,label=false, xticks=(kticks_positioins,pathticks), linestyle = :dashdot)#= ,yrange = [0,2],xrange = [900,2800] =#
end
kticks_positioins
BrillPath[42][:,1]

function GTriang(z, k1, k2)
    numerator = -0.241155 + 0.0000260622 * exp(-16.9014 * z) - 
                0.00186119 * exp(-15.493 * z) + 0.0202915 * exp(-14.0845 * z) - 
                0.104381 * exp(-12.6761 * z) + 0.31285 * exp(-11.2676 * z) - 
                0.576723 * exp(-9.85915 * z) + 0.618961 * exp(-8.4507 * z) - 
                0.262236 * exp(-7.04225 * z) - 0.167707 * exp(-5.6338 * z) + 
                0.217331 * exp(-4.22535 * z) - 0.0999178 * exp(-2.8169 * z) + 
                0.284522 * exp(-1.40845 * z)

    denominator = z + 0.562607 * exp(-16.4384 * z) * (-0.00472295 + 
                 0.0602529 * exp(1.36986 * z) - 0.351124 * exp(2.73973 * z) + 
                 1.22156 * exp(4.10959 * z) - 2.76495 * exp(5.47945 * z) + 
                 4.10544 * exp(6.84932 * z) - 3.67671 * exp(8.21918 * z) + 
                 1.22214 * exp(9.58904 * z) + 1.09802 * exp(10.9589 * z) - 
                 1.37964 * exp(12.3288 * z) + 0.870513 * exp(13.6986 * z) - 
                 1.40078 * exp(15.0685 * z) + exp(16.4384 * z)) * 
                 z * (cos(k1) + 2 * cos(k1 / 2) * cos(sqrt(3) * k2 / 2))

    return -numerator / denominator
end


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







#Calculate the Quantum to Classical Correspondence


function calc_taylorinvmat_fun(corr)
    orderJ,orderω = size(corr)
    @variables x
    taylormat = y -> sum(y[i,1]*x^(i-1) for i =  1:orderJ)
    invtaylormat = inv(taylormat(corr));
    t = Taylor1(Float64,orderJ-1)
    return  substitute(invtaylormat, Dict(x=>t)) 
end

### Compute A 2D Brillouin zone cut: 
N = 40;
kx = (1:N)*2pi/(N);
ky = (1:N)*2pi/(N); #[0.] #for chains
kmat = Tuple.([[1,1/sqrt(3)].*x .+  [1,-1/sqrt(3)].*y for x in kx, y in ky ]);
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);
invstruc = calc_taylorinvmat_fun.(structurefactor);
invcorrs = brillouin_zone_cut_inverse(kmat,invstruc,lattice,center_sites);

(invcorrs[201]/invcorrs[199])

p= plot();
betavec = 0:0.02:10
padetypes = [[4,5],[5,6],[6,6]]
for padetype in padetypes
corrlist = (x->RobustPade.robustpade.(x,padetype[1], padetype[2])).([invcorrs[157],-invcorrs[158],invcorrs[134],invcorrs[159]]);
corrlist = (x->RobustPade.robustpade.(x,padetype[1], padetype[2])).([1/invcorrs[157],-invcorrs[158]/invcorrs[157],invcorrs[134]/invcorrs[157],invcorrs[159]/invcorrs[157]])
for c in corrlist
evaluate(x) = [x(beta) for beta in betavec ] 
corrbetavec = evaluate(c)
plot!(p,betavec,corrbetavec, label = "T = " ,yrange = [-0.1,0.3] )
end
end

sumcoeffs = [sum(abs.(invcorrs[i][:])) for i in 1:length(invcorrs)]
distances = norm.(lattice.sitePositions)
using LaTeXStrings
p = scatter(distances,sumcoeffs, yaxis = :log, title = L"$G_{ij}^{-1} = \sum_n c_n x^n$" , xlabel = "|i-j|" , ylabel = L"\sum_n c_n", label = "Triangular")




sortperm(distances[:])[20]
distances[sortperm(distances[:])][20]

real(1/invcorrs[199])[:]
real(invcorrs[200]/invcorrs[199])[:]
real(invcorrs[223]/invcorrs[199])[:]*10^5
real(invcorrs[201]/invcorrs[199])[:]*10^5
real(invcorrs[135]/invcorrs[199])[:]*10^5