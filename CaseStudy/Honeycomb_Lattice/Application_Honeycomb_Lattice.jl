using JLD2
#using RobustPad
#activates the environment in the folder of the current file

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
lattice,LatGraph,center_sites = getLattice_Ball(L,"honeycomb");
scatter(lattice.sitePositions, xrange = [-10,10],yrange = [-10,10] , aspect_ratio=:equal)

display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.1,fontsize=7,nodeshape=:rect,curves=false))

#2.Compute all correlations in the lattice
#@time Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);
@load "CaseStudy/Honeycomb_Lattice/Correlation_Data_L11.jld2" Correlators




#3. Fourier Transform

### Compute A 2D Brillouin zone cut: 
N = 100
kx = range(-8pi,8pi,length=N)
ky = range(-8pi,8pi,length=N) #[0.] #for chains
kmat = [(y,x) for x in kx, y in ky ]
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
x = -2.5
padetype = [5,6]
evaluate(y) = eval_correlator_LR_continuous_pad(y, x, padetype); #define evaluation function
struc = real.(evaluate.(structurefactor[:,:,1]));
p = Plots.heatmap(kx,ky,struc)
#= ,clims=(0,1.) =#


############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone

#---triangular
path = [(0.01,0.01),(0,2pi/(sqrt(3))),(2pi/(3),2pi/(sqrt(3))),(0.01,0.01)]
pathticks = ["Γ","M","K","Γ"]


#Generate the path 
Nk = 100
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);
#TotalCorrelators=sum(Correlators, dims=2)

###### S(k,w) heatmap
using CairoMakie

x = 2.0
w_vec = collect(0.01:0.0314/2:4.0)
JSkw_mat_total = get_JSkw_mat_finitex("total","pade",x,kvec,w_vec,0.02,1,3,200,false,Correlators,lattice,center_sites)

#plotting
fig = Figure(fontsize=25,resolution = (900,500));
ax=Axis(fig[1,1],limits=(0,Nk+1,0,4),ylabel=L"\omega/J=w",title="Honeycomb: x="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25);
hm=CairoMakie.heatmap!(ax,[k for k in 1:Nk+1],w_vec, JSkw_mat_total,colormap=:viridis,colorrange=(0.001,0.6),highclip=:white);
ax.xticks = (kticks_positioins,pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=40, label = L"J S(k,w)")
resize_to_layout!(fig);
display(fig)

#Add T=0 QMC Data and Linear Spin Wave Theory
disp = [sqrt(1-1/9*abs(exp(1im*k[1]) + exp(1im*(sqrt(3)/2*k[2]-1/2*k[1])) + exp(-1im*(sqrt(3)/2*k[2]+1/2*k[1])))^2)  for k in kvec]
CairoMakie.plot!(ax,[k for k in 1:Nk+1],sqrt(3)*disp, color = :orange, alpha = 0.45,label = L"T=0 \text{LSWT-Dispersion}")

axislegend(ax)
display(fig)

save("HoneycombX"*string(x)*".pdf",fig)
get_pade(m_vec[m_idx],1+Int(floor(max_order/2))-m_idx,1+Int(floor(max_order/2))-m_idx)






#Calculate the Quantum to Classical Correspondence
function calc_taylorinvmat_fun2(corr::Matrix{Matrix{ComplexF64}})
    orderJ,orderω = size(corr[1])
    @variables x::Real
    taylormat = y -> sum(y[i,1]*x^(i-1) for i =  1:orderJ)
    invtaylormat = inv(taylormat.(corr));
    t = Taylor1(ComplexF64,orderJ-1)
    subsfun = y -> substitute(real(y), Dict(x=>t)) + 1im* substitute(imag(y), Dict(x=>t)) 
    taylorinvmat_fun =  subsfun.(invtaylormat)
    return taylorinvmat_fun
end

   #in special cases Num values appear. For those we need to read out their Float value
import Base
function Base.Float64(x::Num)
    return Float64(Symbolics.value(x))
 end     
function calc_taylorinvmat_fun(corr::Matrix{Matrix{ComplexF64}})::Matrix{Taylor1{ComplexF64}}
    orderJ,orderω = size(corr[1])
    @variables x::Real
    taylormat = y -> sum(y[i,1]*x^(i-1) for i =  1:orderJ)
    invtaylormat = inv(taylormat.(corr));
    t = Taylor1(ComplexF64,orderJ-1)
    subsfun = y -> substitute(real(y), Dict(x=>t)) + 1im* substitute(imag(y), Dict(x=>t)) 
    taylorinvmat_fun =  subsfun.(invtaylormat)
    return taylorinvmat_fun
end

### Compute A 2D Brillouin zone cut: 
N = 30;
kx = (1:N)*2pi/(N);
ky = (1:N)*2pi/(N); #[0.] #for chains
kmat = Tuple.([[1,1/sqrt(3)].*x .+ [1,-1/sqrt(3)].*y for x in kx, y in ky ]);
structurefactor =  brillouin_zone_cut_Matrix(kmat,Correlators,lattice,center_sites);
invstruc = Array{Matrix{Taylor1{ComplexF64}}}(undef,N,N)
Threads.@threads for i=1:N
    for j = 1:N
        invstruc[i,j] = calc_taylorinvmat_fun(structurefactor[i,j])
    end
end
invcorrs = inverse_fourier_transform(kmat,invstruc,lattice,center_sites);


structurefactor =  brillouin_zone_cut_Matrix(kmat,Correlators,lattice,center_sites);
invcorrstest = inverse_fourier_transform(kmat,structurefactor,lattice,center_sites);

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



sumcoeffs = [sum(abs.(invcorrs[i,b][:])) for i in 1:length(lattice), b in 1:2]
distances = [norm(lattice.sitePositions[i] .- lattice.sitePositions[center_sites[b]]) for i in 1:length(lattice) , b in 1:2]
pHoneycomb = scatter(pkagome,distances[:],sumcoeffs[:], yaxis = :log, title = L"$G_{ij}^{-1} = \sum_n c_n x^n$" , xlabel = "|i-j|" , ylabel = L"\sum_n |c_n|", label = "Honeycomb")



sortperm(distances[:,1])[14]
distances[sortperm(distances[:,1]),1][14]


real(1/invcorrs[108,1])[:]
real(invcorrs[109,1]/invcorrs[108,1])[:]
real(invcorrs[84,1]/invcorrs[108,1])[:]*10^5
real(invcorrs[83,1]/invcorrs[108,1])[:]*10^5
real(invcorrs[63,1]/invcorrs[108,1])[:]*10^5
