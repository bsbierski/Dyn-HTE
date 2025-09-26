using JLD2,DelimitedFiles,Polynomials
using DynHTE
include("../plotConventions.jl")

### load graph evaluations and prepare lattice  
L = 12
n_max = 1*L
spin_length = 1/2
hte_graphs = load_dyn_hte_graphs(spin_length,L);
hte_lattice = getLattice(L,"square");

### plot lattice if desired
#display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

### compute all correlations in the lattice by graph embedding (or load them)
fileName_c = "examples/Square_Lattice/Square_Lattice_"*create_spin_string(spin_length)*"_c_iipDyn_nmax"*string(n_max)*"_L"*string(L)*".jld2"
if isfile(fileName_c)
    c_iipDyn_mat = load_object(fileName_c)
else
    c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs);
    save_object(fileName_c,c_iipDyn_mat)
end


#########################################################################################
###### DSF linecut for single k #########################################################
#########################################################################################

### pick one of the three high-symmetry k-points
k,k_label = (π/2,π/2), "(piO2,piO2)" 
#k,k_label = (1.0*π,0.0), "(pi,0)"

### compute moments
c_kDyn = get_c_k(k,c_iipDyn_mat,hte_lattice)
m_vec = get_moments_from_c_kDyn(c_kDyn)

### define x vectors for plotting moments, for the bare series and the few x0 where 
x_vec_bare = collect(0:0.025:1.4)
x_vec = collect(0.0:0.1:2.5)
x0_vec = [1.0,1.5,2.0] 
m0_vec = [Float64[] for _ in x0_vec]

### re-sum moments with Padé approximant of u-series
# maximal order for moment
r_max = 3  

# f-factor for transform to u-series
f=0.7
ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)
u_vec = tanh.(f .* x_vec)
u0_vec = tanh.(f .* x0_vec)

### plot of moments over x
plt_m = Plots.plot([0],[0],xlims=(0,x_vec[end]),label="",xlabel=L"x=J/T",ylabel=L"x \cdot m_{\mathbf{k},2r}(x) \, / \, m_{\mathbf{k},2r}(0)",legend=:topleft);
Plots.plot!(plt_m,-x_vec,x_vec,color=:grey,label="x bare");
Plots.plot!(plt_m,-x_vec,x_vec,color=:grey,linestyle=linestyle_vec[2],label="u Padé [7-r,6-r]");
Plots.plot!(plt_m,-x_vec,x_vec,color=:grey,linestyle=linestyle_vec[3],label="u Padé [6-r,5-r]");
Plots.annotate!(plt_m,1.5,1,L"$\mathbf{k}=$"*string(k_label)*",  f="*string(f),7);

Plots.plot!(plt_m,title="SquareLattice AFM S=1/2, k="*k_label*" (f=$f)")
for r in 0:r_max
    #xm_norm_r = m_vec[1+r]
    xm_norm_r = coeffs(Polynomial([0,1],:x) * (m_vec[1+r]/m_vec[1+r](0)))
    println()
    println("r=$r")
    p_u = Polynomial(ufromx_mat[1:n_max+2-2*r,1:n_max+2-2*r]*xm_norm_r)
    
    Plots.plot!(plt_m,x_vec_bare,Polynomial(xm_norm_r).(x_vec_bare),color=color_vec[r+1],linewidth=0.4,label="r=$r",alpha=0.7) 
    #plot!(plt_m,x_vec,p_u.(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[1],label="u-poly r=$r")
    
    Plots.plot!(plt_m,x_vec,get_pade(p_u,7-r,6-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="",alpha=0.7)
    Plots.plot!(plt_m,x_vec,get_pade(p_u,6-r,5-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="",alpha=0.7)

    ### extract moments at x0_vec 
    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]
        u0 = u0_vec[x0_pos]
        push!(m0_vec[x0_pos],m_vec[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
    end
end

### For x0 ∈ x0_vec: Plot the continued fraction parameters δ_r and the DSF
w_vec = collect(0.0:0.025:5)
plt_JS = Plots.plot([0,0],[0,0.01],color=:grey,label="Dyn-HTE",xlims=(0,w_vec[end]),xlabel=L"\omega/J=w",ylabel=L"J \, S(\mathbf{k}="*k_label*L",\omega)",legend=:topleft);#,title="SquareLat AFM S=1/2: "*"JS("*k_label*",w)")
plt_δ =Plots.plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_{\mathbf{k},r}",legend=:topleft);

for x0_pos in eachindex(x0_vec)
    x0 = x0_vec[x0_pos]

    ### plot Dyn-HTE
    δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
    Plots.scatter!(plt_δ,r_vec,δ_vec,color=thermalCol4_vec[x0_pos],label="x=$x0")
    Plots.plot!(plt_JS,w_vec,[JSwithTerminator(δ_vec,x0,w,get_extrapolation_params(δ_vec,3,3,true)) for w in w_vec],color=thermalCol4_vec[x0_pos],label="")

    ### plot experimental data from DallaPiazza for comparison, adjust factor for relative scaling if needed
    fileName = "examples/Square_Lattice/DallaPiazza_exp_S"*k_label*".csv"
    factor = 1200
    if isfile(fileName)
        data = readdlm(fileName,',',Float64)
        Plots.scatter!(plt_JS,data[:,1],data[:,2]/factor,color=:blue,marker=:cross,markersize=3.0,label="")#,label="x=$x0"
    end
end

### combine the plots
xPlots,yPlots=1,3
plt_final = Plots.plot(plt_m,plt_δ,plt_JS,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
display(plt_final)


#########################################################################################
###### DSF for k-path ###################################################################
#########################################################################################

### define plot paramters
w_vec = collect(0.0:0.025:3.5)
r_max = 3                
f=0.7
ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)
x0 = 2.0
u0 = tanh.(f .* x0)

### define and generate k-path 
path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
pathticks = ["(π/2,π/2)","(π,0)","(π,π)","(π/2,π/2)","(0,0)","(π,0)"]

Nk = 75
k_vec,kticks_positioins = create_brillouin_zone_path(path, Nk)
JSkw_mat = zeros(Nk+1,length(w_vec))

### fill JSkw_mat
for k_pos in eachindex(k_vec)
    k = k_vec[k_pos]

    c_kDyn = get_c_k(k,c_iipDyn_mat,hte_lattice)
    m_vec = get_moments_from_c_kDyn(c_kDyn)
    m0 = Float64[]

    for r in 0:r_max
        xm_norm_r = coeffs(Polynomial([0,1],:x) * (m_vec[1+r]/m_vec[1+r](0)))
        p_u = Polynomial(ufromx_mat[1:n_max+2-2*r,1:n_max+2-2*r]*xm_norm_r)
        push!(m0,m_vec[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
    end

    δ_vec,r_vec = fromMomentsToδ(m0)
    JSkw_mat[k_pos,:] = [JSwithTerminator(δ_vec,x0,w,get_extrapolation_params(δ_vec,r_max,r_max,true)) for w in w_vec]#[JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec]
end

### plot JS(k,ω)
using CairoMakie

fig = CairoMakie.Figure(fontsize=8,size=(aps_width,0.6*aps_width));
ax=Axis(fig[1,1],xlabel=L"\mathbf{k}",ylabel=L"\omega/J=w",xlabelsize=8,ylabelsize=8);
hm=CairoMakie.heatmap!(ax,collect(0:Nk)/(Nk),w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,0.4),highclip=:white);
ax.xticks = ((kticks_positioins .- 1)/(Nk),pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=8, label = L"J S(\mathbf{k},\omega)")
CairoMakie.text!(ax,"x=J/T=$x0",position=[(0.05,0.5)],color=:white)
CairoMakie.text!(ax,"f=$f",position=[(0.05,0.2)],color=:white)

#Add experimental magnon data from DallaPiazza at T/J~ 0.09
data = readdlm("examples/Square_Lattice/DallaPiazza_exp_MagnonEnergy.csv",',',Float64)
CairoMakie.plot!(ax,data[:,1],data[:,2],color=:purple,marker=:cross,markersize=5.0,label="")

resize_to_layout!(fig);
display(fig)