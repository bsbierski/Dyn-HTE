using JLD2
#using RobustPad

include("../../Embedding.jl")
include("../../LatticeGraphs.jl")
include("../../ConvenienceFunctions.jl") 
#specify max order
max_order = 12

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
lattice,LatGraph,center_sites = getLattice_Ball(L,"square");
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.12,fontsize=4,nodeshape=:rect,curves=false))

#2.Compute or load all correlations in the lattice
#@time Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);
@load "CaseStudy/Square_Lattice/Correlation_Data_L12.jld2" Correlators

#check the susceptibility with 10.1103/PhysRevB.53.14228
(brillouin_zone_cut([(0.0,0.0) (0.0,0.0) ;(0.0,0.0)  (0.0,0.0)],Correlators,lattice,center_sites)[1]*4)[:,1].*[factorial(n)*4^n for n in 0:max_order]

Correlators

#3. Fourier Transform

### Compute A 2D Brillouin zone cut: 
N = 30
kx = [0.,0]#(1:N)*2pi/(N)
ky = [0.,0]#(1:N)*2pi/(N) #[0.] #for chains
kmat = [(y,x) for x in kx, y in ky ]
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites)
corrs=brillouin_zone_cut_inverse(kmat,structurefactor,lattice,center_sites);


structurefactor[1][:,1]
corrs[1]
Float64.(Correlators[1])


structurefactor[1,1]

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
x = -1.2
padetype = [5,6]
evaluate(y) = eval_correlator_LR_continuous_pad(y, x, padetype); #define evaluation function
struc = real.(evaluate.(structurefactor));
p = Plots.heatmap(kx,ky,struc,clims=(0,1.))



### Calculate the full S(ω,k)
x = 2.01
N = 100
kx = range(-0+0.001,2pi+0.001,length=N)
ky = range(-0+0.001,2pi+0.001,length=N) 
w_vec = collect(-5:10/100:5.0)
struc = zeros(N,N,length(w_vec));
Threads.@threads for i = 1:N
for  j = 1:N
    struc[i,j,:] .=  get_JSkw_mat_finitex("total","pade",x,[[kx[i],kx[j]]],w_vec,0.02,1,3,200,false,Correlators,lattice,center_sites)[1,:]
end
end
#calculate the integral identity
4*sum(struc)*((maximum(w_vec)-minimum(w_vec))/length(w_vec))*(4*pi^2/N^2)/(4*pi^2)

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
w_ind = 75
Plots.heatmap(kx,ky,struc[:,:,w_ind], title = L"$S(\omega,k)$ at $\omega$ ="*string(w_vec[w_ind]))




####Plot the susceptibility along the high symmetry path for various temperatures 
Tvec = [2.0,1.0,0.5,0.375]
padetype = [4,4]

p= Plots.plot();
for T in Tvec
evaluate(x) = eval_correlator_LR_continuous_pad(x,-1/T,padetype);
cut = 1/T* real.(evaluate.(BrillPath));
Plots.plot!(p,1:length(kvec),cut, label = "T = $T", xticks=(kticks_positioins,pathticks)#= ,yrange = [0,2],xrange = [900,2800] =#)
end
p




############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone
#---square
path = [(0.001,0.001),(pi,0),(pi,pi),(0.001,0.001)]
pathticks = ["Γ","X","M","Γ"]

#---square path Sherman
path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
pathticks = ["K","X","M","K","Γ","X"]


#Generate the path 
Nk = 100
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);
kticks_positioins
BrillPath[30][:,1]
###### S(k,w) heatmap
using CairoMakie

x = 2.5
w_vec = collect(0.01:0.0314/2:5.0)
JSkw_mat = get_JSkw_mat_finitex("total","pade",x,kvec,w_vec,0.00,1,2,200,false,Correlators,lattice,center_sites)


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

inv(10)




calc_taylorinvmat_fun(corr)


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
N = 30;
kx = (1:N)*2pi/(N);
ky = (1:N)*2pi/(N); #[0.] #for chains
kmat = [(y,x) for x in kx, y in ky ];
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);
invstruc = calc_taylorinvmat_fun.(structurefactor);
invcorrs = brillouin_zone_cut_inverse(kmat,invstruc,lattice,center_sites);


#evaluate first then invert 
eval(x) = eval_correlator_LR_continuous_pad_tanh(x,-1/0.2,[6,6],2);

invstruc = 1 ./ eval.(structurefactor)
invcorrs = brillouin_zone_cut_inverse(kmat,invstruc ,lattice,center_sites);

fnew = []
gnew = []
ϵ1new = []
ϵ2new = []
betas = 1:5
for beta = betas
    eval(x) = eval_correlator_LR_continuous_pad_tanh(x,-beta,[12,0],0.8);
    invstruc = 1 ./ eval.(structurefactor)
    invcorrs = brillouin_zone_cut_inverse(kmat,invstruc ,lattice,center_sites);
    append!(fnew,real(1/invcorrs[157]))
    append!(gnew,real(invcorrs[181,1]/invcorrs[157,1]))
    append!(ϵ1new,real(invcorrs[132,1]/invcorrs[157,1]))
    append!(ϵ2new,real(invcorrs[111,1]/invcorrs[157,1]))
end
p= plot(betas,fnew, ylims = [0,0.3]);
plot!(p,betas,gnew);
plot!(p,betas,ϵ1new);
plot!(p,betas,ϵ2new)

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
psquare=scatter(p,distances,sumcoeffs, yaxis = :log, title = L"$G_{ij}^{-1}(i\omega = 0) = \sum_n c_{ij,n} x^n$" , xlabel = "|i-j|" , ylabel = L"\sum_n |c_{ij,n}|", label = "Square")




sortperm(distances[:])[14]
distances[sortperm(distances[:])][14]

real(1/invcorrs[157])[:]
real(invcorrs[181,1]/invcorrs[157,1])[:]
real(invcorrs[132,1]/invcorrs[157,1])[:]*10^5
real(invcorrs[111,1]/invcorrs[157,1])[:]*10^5
real(invcorrs[110,1]/invcorrs[157,1])[:]*10^5



# QMC data for QCC Square lattice
using HDF5
fid = h5open("CaseStudy/Square_Lattice/job_BSb_beta1.out.h5", "r")
res=fid["simulation"]["results"]["DensDens_CorrFun_w0"]["mean"]["value"][:]
L = Int(length(res)^(1/2))
positions = [(n -L*(floor(n/(L/2))), m - L*floor(m/(L/2)) ) for  n=0:(L-1) , m=0:(L-1) ]
pos = reshape(positions,L^2)
resmat= reshape(res,L,L)
ressym = 1/2*(resmat+transpose(resmat))[:]
N = Int(40);
kx = (1:N)*2pi/(N);
ky = (1:N)*2pi/(N); #[0.] #for chains
structurefactor1 =  brillouin_zone_cut([(y,x) for x in kx, y in ky ],ressym,positions);
structurefactor2 =  brillouin_zone_cut([(y,-x) for x in kx, y in ky ],ressym,positions);
structurefactor3 =  brillouin_zone_cut([(-y,x) for x in kx, y in ky ],ressym,positions);
structurefactor4 =  brillouin_zone_cut([(-y,-x) for x in kx, y in ky ],ressym,positions);
structurefactor = 1/4*(structurefactor1+structurefactor2+structurefactor3+structurefactor4)
minimum(structurefactor)


Plots.heatmap(structurefactor)

invtest
res

invstruc = 1 ./structurefactor
invcorrs = brillouin_zone_cut_inverse(kmat,invstruc ,lattice,center_sites);
Plots.heatmap(structurefactor)



f = []
g = []
ϵ1 = []
ϵ2 = []
betas = 1:5
for beta = betas
    fid = h5open("CaseStudy/Square_Lattice/job_BSb_beta$beta.out.h5", "r")
    res=fid["simulation"]["results"]["DensDens_CorrFun_w0"]["mean"]["value"][:]
    L = Int(length(res)^(1/2))
    positions = [(n -L*(floor(n/(L/2))), m - L*floor(m/(L/2)) ) for  n=0:(L-1) , m=0:(L-1) ]
    pos = reshape(positions,L^2)
    resmat= reshape(res,L,L)
    ressym = 1/2*(resmat+transpose(resmat))[:]
    N = L;
    kx = (1:N)*2pi/(N);
    ky = (1:N)*2pi/(N); #[0.] #for chains
    structurefactor1 =  brillouin_zone_cut([(y,x) for x in kx, y in ky ],ressym,positions);
    structurefactor2 =  brillouin_zone_cut([(y,-x) for x in kx, y in ky ],ressym,positions);
    structurefactor3 =  brillouin_zone_cut([(-y,x) for x in kx, y in ky ],ressym,positions);
    structurefactor4 =  brillouin_zone_cut([(-y,-x) for x in kx, y in ky ],ressym,positions);
    structurefactor = 1/4*(structurefactor1+structurefactor2+structurefactor3+structurefactor4)
    invstruc = 1 ./structurefactor
    invcorrs = brillouin_zone_cut_inverse2(kmat,invstruc,positions);
    append!(f,real(1/invcorrs[1]))
    append!(g,real(invcorrs[2]/invcorrs[1]))
    append!(ϵ1,real(invcorrs[3]/invcorrs[1]))
    append!(ϵ2,real(invcorrs[L+2]/invcorrs[1]))
end
p= plot(betas,f , ylims = [0,0.3]);
plot!(p,betas,g);
plot!(p,betas,ϵ1);
plot!(p,betas,ϵ2)
scatter!(p,betas,fnew)


structurefactor =  brillouin_zone_cut(kmat,res[sortperm(norm.(pos))][1:9],pos[sortperm(norm.(pos))][1:9]);

scatter(norm.(pos[sortperm(norm.(pos))]),abs.(res[sortperm(norm.(pos))]) , xlims = [0,10] )

norm1(x) = norm(x,1) 


scatter(positions[norm1.(positions) .< 10], xrange = [-10,10],yrange = [-10,10] , aspect_ratio=:equal)


position
norm.(position) == transpose(norm.(position))
resmat= reshape(res,L,L)

resmatsym = (1/2*(resmat+transpose(resmat)))[:]

function brillouin_zone_cut(kmat::Union{Matrix{Tuple{Float64,Float64}},Matrix{Tuple{Float64,Float64,Float64}}},Correlators::Vector{Float64},positions)
    """computes the fourier transform along a 2D cut through the 2D or 3D k-space
        given the Correlation Matrix computet from compute_lattice_correlations """
    (nx,ny) = size(kmat)

    structurefactor = Array{Float64}(undef, nx,ny);
    Threads.@threads for i in 1:nx
        for j in 1:ny
            z = 0
                for k in 1:length(positions)
                    if norm(positions[k],1) < 30
                    z += cos(dot(kmat[i,j], positions[k])) *  Correlators[k]
                    end
                 end
            structurefactor[i,j] = abs(z) #= / (length(lattice) * length(lattice.unitcell.basis)) =#
        end
    end
    return structurefactor
end


function brillouin_zone_cut_inverse2(kmat::Union{Matrix{Tuple{Float64,Float64}},Matrix{Tuple{Float64,Float64,Float64}}},Correlators,positions)
    """computes the fourier transform along a 2D cut through the 2D or 3D k-space
        given the Correlation Matrix computet from compute_lattice_correlations """
    (nx,ny) = size(kmat)

    structurefactor = Array{Float64}(undef, length(positions));
    Threads.@threads for k in 1:length(positions)
        z = 0

     for i in 1:nx
        for j in 1:ny
                    z += cos(dot(kmat[i,j], positions[k])) *  Correlators[i,j]
            end 
        end
    structurefactor[k] = z/(nx*ny) #= / (length(lattice) * length(lattice.unitcell.basis)) =#
end
    return structurefactor
end


norm((1,1))