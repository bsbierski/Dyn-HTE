using JLD2
#using RobustPad
#activates the environment in the folder of the current file

include("../../LatticeGraphs.jl")
include("../../Embedding.jl")
include("../../ConvenienceFunctions.jl") 
#specify max order


L = 0
spin_length = 1/2
hte_graphs = load_dyn_hte_graphs(spin_length,L);
hte_lattice = getLattice(L,"kagome");
display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))
@time c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs);

hte_lattice.lattice.sitePositions
scatter(hte_lattice.lattice.sitePositions)
scatter!([(0.5, 0.8660254037844387)])
#test uniform susceptibility with /10.1103/PhysRevB.89.014415 
(brillouin_zone_cut([(0.0,0.0) (0.0,0.0) ;(0.0,0.0)  (0.0,0.0)],Correlators,lattice,center_sites)[1]*4/3)[:,1].*[factorial(n)*4^n for n in 0:max_order]


#3. Fourier Transform
### Compute A 2D Brillouin zone cut: }
N = 40
kx = range(-2pi,2pi,length=N)
ky = range(-2pi,2pi,length=N) #[0.] #for chains
kmat = [(y,x) for x in kx, y in ky ]
c_kDyn =  get_c_kDyn(kmat,c_iipDyn_mat,hte_lattice);
### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
x = -2
padetype = [2,2]
evaluate(y) = eval_correlator_LR_continuous_pad(y, x, padetype); #define evaluation function
struc = ( evaluate.(c_kDyn));
p = Plots.heatmap(kx,ky,struc, clims=(0,0.3))


### Calculate the full S(ω,k)
x = 2.0
N = 20
kx = range(-0+0.001,2pi+0.001,length=N)
ky = range(-0+0.001,2pi+0.001,length=N) 
w_vec = collect(-5:10/10:5.0)
kmat = [[1,1/sqrt(3)].*x .+  [1,-1/sqrt(3)].*y for x in kx, y in ky ]
Swk = zeros(N,N,length(w_vec));
Threads.@threads for i = 1:N
for  j = 1:N
    Swk[i,j,:] .=  get_JSkw_mat_finitex("total","pade",x,[kmat[i,j]],w_vec,0.02,1,3,200,false,Correlators,lattice,center_sites)[1,:]
end
end

##Check sum rule
4/3*sum(Swk)*((maximum(w_vec)-minimum(w_vec))/length(w_vec))*(4*pi^2/N^2)/(4*pi^2)

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
w_ind = 4
Plots.heatmap(kx,ky,Swk[:,:,w_ind], title = L"$S(\omega,k)$ at $\omega$ ="*string(w_vec[w_ind]))



############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone

#---triangular
path = [(0.1,0.1),(0,2pi/(sqrt(3))),(2pi/(3),2pi/(sqrt(3))),(0.1,0.1)]
pathticks = ["Γ","M","K","Γ"]


#Generate the path 
Nk = 100
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);
 
kticks_positioins
JSkw_mat_total = get_JSkw_mat_finitex("total","pade",x,[kvec[37]],w_vec,0.02,1,3,200,false,Correlators,lattice,center_sites)
###### S(k,w) heatmap
using CairoMakie

x = 1.5
w_vec = collect(0.01:0.0314/2:4.0)
JSkw_mat_total = get_JSkw_mat_finitex("total","pade",x,kvec,w_vec,0.02,1,3,200,false,Correlators,lattice,center_sites)

#plotting
fig = Figure(fontsize=25,resolution = (900,500));
ax=Axis(fig[1,1],limits=(0,Nk+1,0,4),ylabel=L"\omega/J=w",title="Kagome: x="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25);
hm=CairoMakie.heatmap!(ax,[k for k in 1:Nk+1],w_vec, JSkw_mat_total,colormap=:viridis,colorrange=(0.001,1.0),highclip=:white);
ax.xticks = (kticks_positioins,pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=40, label = L"J S(k,w)")
resize_to_layout!(fig);
display(fig)

save("./KagomeX"*string(x)*".pdf",fig)

taylor(::Polynomial, ::Any...; kwargs...)


function c_n(n, r)
    if n == 1
        return (1/3) * r
    elseif n == 2
        return -(4/9) * r^2
    elseif n == 3
        return (1/9) * r^2 * (-1 + 4r)
    elseif n == 4
        return -(4/405) * r^2 * (3 - 28r + 37r^2)
    elseif n == 5
        return (1/4860) * r^2 * (-45 + 702r - 1892r^2 + 1328r^3)
    elseif n == 6
        return -(1/510300) * r^2 * (1728 - 35946r + 164289r^2 - 207896r^3 + 102576r^4)
    elseif n == 7
        return (1/6123600) * r^2 * (-8694 + 218916r - 1401381r^2 + 2888772r^3 - 2251248r^4 + 909184r^5)
    elseif n == 8
        return -(1/22963500) * r^2 * (15390 - 446256r + 3538764r^2 - 10535337r^3 + 12202552r^4 - 7318640r^5 + 2416640r^6)
    elseif n == 9
        return (1/7715736000) * r^2 * (-2710665 + 87954822r - 807482331r^2 + 3091042674r^3 - 5118502560r^4 + 4009481184r^5 - 2113197952r^6 + 518354176r^7)
    else
        return 0  # Undefined for n > 9 in given series
    end
end

(brillouin_zone_cut([(0.0,0.0) (0.0,0.0) ;(0.0,0.0)  (0.0,0.0)],Correlators,lattice,center_sites)[1]/3)[:,1]


[c_n(n, 1/2*(1/2+1)) for n = 1:8]





#3. Fourier Transform
### Compute A 2D Brillouin zone cut by calculating the pade for the inverse matrix}
function calc_taylorinvmat_fun(corr)
    orderJ,orderω = size(corr[1])
    @variables x
    taylormat = y -> sum(y[i,1]*x^(i-1) for i =  1:orderJ)
    invtaylormat = inv(taylormat.(corr));
    t = Taylor1(Float64,orderJ)
    return taylorinvmat_fun = substitute(invtaylormat, Dict(x=>t)) 
end

function eval(taylorinvmat_fun,X, padetype)
    fun = f -> RobustPade.robustpade.(f,padetype[1], padetype[2])
    mat = (fun.(taylorinvmat_fun))
    return sum(inv(map(f -> f(X), mat)))
end

N = 
kx = range(-2pi,2pi,length=N)
ky = range(-2pi,2pi,length=N) #[0.] #for chains
kmat = [(y,x) for x in kx, y in ky ]
structurefactor =  brillouin_zone_cut_Matrix(kmat,Correlators,lattice,center_sites);
### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
taylor_vec = Vector{Matrix{Taylor1{Float64}}}(undef,N^2);
for i = 1:N^2
    taylor_vec[i] = calc_taylorinvmat_fun(structurefactor[i])
end

X = -3
padetype = [6,5]
evaluate(y) = eval(y, X, padetype); #define evaluation function
res_vec = Vector{Float64}(undef,N^2);
Threads.@threads for i = 1:N^2
    res_vec[i] = evaluate(taylor_vec[i])
end
struc = reshape(res_vec,(N,N))
p = Plots.heatmap(kx,ky,struc, clims=(0,1.))

#Calculate the Quantum to Classical Correspondence

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

kmat
### Compute A 2D Brillouin zone cut: 
N = 20;
kx = (1:N)*2pi/(N);
ky = (1:N)*2pi/(N); #[0.] #for chains
kmat = Tuple.([[1,1/sqrt(3)].*x .+ [1,-1/sqrt(3)].*y for x in kx, y in ky ]);
structurefactor =  get_c_kDyn_subl(kmat,c_iipDyn_mat,hte_lattice);
invstruc = Array{Matrix{Taylor1{ComplexF64}}}(undef,N,N)
Threads.@threads for i=1:N
    for j = 1:N
        invstruc[i,j] = calc_taylorinvmat_fun(structurefactor[i,j])

    end
end
invcorrs = inverse_fourier_transform(kmat,invstruc,lattice,center_sites);


structurefactor =  get_c_kDyn_subl(kmat,c_iipDyn_mat,hte_lattice);
invcorrstest = inverse_fourier_transform_subl(kmat,structurefactor,hte_lattice);

hte_lattice.basis_positions
invcorrstest[132]
Float64.(c_iipDyn_mat[132])



typeof(structurefactor)
structurefactor[1]
(1/invcorrs[190,1])[:]
Float64.(Correlators[193])
invcorrstest[193]

(invcorrs[201]/invcorrs[199])

center_sites

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


center_sites
real(1/invcorrs[190,1])[:]
real(invcorrs[190,2]/invcorrs[190,1])[:]
real(invcorrs[193,1]/invcorrs[190,1])[:]*10^5



sumcoeffs = [sum(abs.(invcorrs[i,b][:])) for i in 1:length(lattice), b in 1:3]
distances = [norm(lattice.sitePositions[i] .- lattice.sitePositions[center_sites[b]]) for i in 1:length(lattice) , b in 1:3]
using LaTeXStrings
pkagome = scatter(distances[:],sumcoeffs[:], yaxis = :log, title = L"$G_{ij}^{-1} = \sum_n c_n x^n$" , xlabel = "|i-j|" , ylabel = L"\sum_n c_n", label = "Kagome")


sortperm(distances[:,1])[16]
distances[sortperm(distances[:,1])][16]

real(1/invcorrs[190,1])[:]
real(invcorrs[191,1]/invcorrs[190,1])[:]
real(invcorrs[228,1]/invcorrs[190,1])[:]*10^5
real(invcorrs[193,1]/invcorrs[190,1])[:]*10^5
real(invcorrs[194,1]/invcorrs[190,1])[:]*10^5


sort()

sumcoeffs[:]