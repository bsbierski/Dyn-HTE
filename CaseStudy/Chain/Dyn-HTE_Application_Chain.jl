######### Dyn-HTE for chain #########
using Polynomials, HDF5, Measurements


path_DynHTSE="C:/Users/ruben/Documents/GitHub/Projects/Master/Dyn-HTE/Dyn-HTE/"
include(path_DynHTSE*"LatticeGraphs.jl")
include(path_DynHTSE*"ConvenienceFunctions.jl")
include(path_DynHTSE*"plotConventions.jl") 
include(path_DynHTSE*"Embedding.jl")


#DEFINE THE SYSTEM
n_max = 12                   #max order in perturbation theory (currently 12 are available)
spin_length = 1/2            #spin length (currently only S=1/2,1 are available)
#load graphs
hte_graphs = load_dyn_hte_graphs(spin_length,n_max)
#define the lattice geometry
hte_lattice = getLattice(n_max,"chain")
#calculate the correlation function for all site combinations
c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs)




#### moments and δs of the continued fraction expansion

#we will calculate the first four moments at fixed k

k= 0.2*pi            #define fixed k 
x_vec = 0:0.01:4.5   #define temperature range of interest
#Fourier transform the correlation functions at k
c_kDyn_mat = get_c_k([(k,0.0)],c_iipDyn_mat,hte_lattice)[1]
#calculate the moments 
m_vec = get_moments_from_c_kDyn(c_kDyn_mat)
#rescale moments 
m_vec_times_x = [m_vec[i]*Polynomial([0,1]) for i=1:length(m_vec)]

##EXTRAPLOATION OF MOMENTS

#basic pade 
m_vec_extrapolated_pade = []
for m_idx=1:length(m_vec)-2
    push!(m_vec_extrapolated_pade   , extrapolate_series(m_vec[m_idx],"pade",(7-m_idx,7-m_idx)))
end

#pade in u = tanh(f*x) (2 different versions)
f= 0.48   #define the f value (f=0.48 is very fine tuned to give good results)
m_vec_times_x_extrapolated_u_pade1 = []
m_vec_times_x_extrapolated_u_pade2 = []
for m_idx=1:length(m_vec)-2
    push!(m_vec_times_x_extrapolated_u_pade1, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(8-m_idx,7-m_idx,f)))
    push!(m_vec_times_x_extrapolated_u_pade2, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(7-m_idx,8-m_idx,f)))
end

#plot the moments 
plt_m = Plots.plot([0],[0],label="",xlabel=L"x",ylabel=L"x \cdot m_{\mathbf{k},r}(x)/m_{\mathbf{k},r}(0)",title="moments at k="*string(k/pi)*"π",legend=:topleft,xlim=(-0.2,4.5),size=(1.5*aps_width,aps_width))
Plots.plot!(plt_m,[0],[0],label="x bare",color = "grey",linestyle = linestyle_vec[1],linewidth=0.4)
Plots.plot!(plt_m,[0],[0],label="u Padé [7-r,6-r]",color = "grey",linestyle = linestyle_vec[2],alpha =0.5)
Plots.plot!(plt_m,[0],[0],label="u Padé [6-r,7-r]",color = "grey",linestyle = linestyle_vec[3])
Plots.plot!(plt_m,[0],[0],label="Padé [6-r,6-r]",color = "grey",linestyle = linestyle_vec[4])
for i=1:4
    Plots.plot!(plt_m,[0],[0],label="r="*string(i-1),color = color_vec[i])
    Plots.plot!(plt_m,x_vec[1:180],m_vec_times_x[i].(x_vec[1:180])./(m_vec[i](0)),label = "",alpha= 0.7,color = color_vec[i],linestyle = linestyle_vec[1],linewidth=0.5)
    Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade1[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 0.5,color = color_vec[i],linestyle = linestyle_vec[2],linewidth=0.5)
    Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade2[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[3],linewidth=0.5)
    Plots.plot!(plt_m,x_vec,x_vec.*m_vec_extrapolated_pade[i].(x_vec)./(m_vec_extrapolated_pade[i](0.0001)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[4],linewidth=0.5)
end

display(plt_m)

#now calculate the δs of the continued fraction expansion 
#we will use the u_pade extrapolation

x_vec = [0.5,1.0,2.0,4.0]  #define the x values at which we calculate the deltas

deltas_at_x =[]
for x in x_vec
    δ_vec,r_vec = fromMomentsToδ([m(tanh(f*x))/x for m in m_vec_times_x_extrapolated_u_pade1])
    push!(deltas_at_x,δ_vec)
end

#plots the δs 
plt_δ = Plots.plot([0],[0],title="δs using Padé with substitution on moments",label="",xlabel=L"r",ylabel=L"\delta_r",legend=:topleft,size=(1.5*aps_width,aps_width))
for (i,x) in enumerate(x_vec)
    Plots.scatter!(plt_δ,0:3,deltas_at_x[i][1:4],label = "x="*string(x),color = thermalCol4_vec[i],markersize=7)
end

display(plt_δ)


#one can now go on and extrapolate the deltas to infinity via
deltas_at_x_ext = []
for x_idx = 1:length(x_vec)
    δ_vec_ext = extrapolate_δvec(deltas_at_x[x_idx],3,3,1000,false)
    push!(deltas_at_x_ext,δ_vec_ext)
end

#add the extrapolated deltas to the plot 
for (i,x) in enumerate(x_vec)
    Plots.scatter!(plt_δ,4:10,deltas_at_x_ext[i][5:11],label = "",title="δs extrapolated",color = thermalCol4_vec[i],markersize=4,alpha=0.5)
end

display(plt_δ)

#we can use the deltas to calculate the spin structure factor (still at fixed k)
w_vec = -3:0.01:3  #define w
x_idx = 3          #choose temperature index (relativ to x_vec = [0.5,1.0,2.0,4.0]) 
η = 0.01           #the imaginary part after analytic continuation (η ->0)
#now calculate the DSF 
DSF = [JS(deltas_at_x_ext[x_idx] ,x_vec[x_idx],w,η) for w in w_vec]

plt_dsf = Plots.plot(w_vec,DSF,label="",title="DSF at k="*string(k/pi)*"π and x="*string(x_vec[x_idx]),xlabel=L"\omega/J",ylabel=L"JS(k,\omega)",size = (1.5*aps_width,aps_width))

display(plt_dsf)


###SPIN STRUCTURE FACTOR HEATMAPS
using CairoMakie

x = 0.0                #define temperature (x=J/T)
k_step_size = 1/41     #define k step size (in 1/π)
w_step_size = 0.025    #define ω step size (in 1/J)
#define k and ω vectors 
k_vec = vcat(vcat((0.0001,0.0),[(k*pi,0.0) for k in 0:k_step_size:2][2:end-1] ),(1.999*pi,0.0))#[(k,0.0) for k in 0.01:0.0039*2.4:(2*π-0.01)]
w_vec = collect(-3:w_step_size:3)

#calculate the spin structure factor for the given k and ω 
JSkw_mat = get_JSkw_mat("u_pade",x,k_vec,w_vec,c_iipDyn_mat,hte_lattice)

#plot the result
fig = Figure(size=(400,400),fontsize=20)
ax=Axis(fig[1,1],limits=(0,2,-3,3),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title="x="*string(x),titlesize=20,xlabelsize=20,ylabelsize=20)
hm=CairoMakie.heatmap!(ax,[k[1]/π for k in k_vec],w_vec,JSkw_mat,colormap=:viridis,colorrange=(0.0,0.45),highclip=:white)
display(fig)

