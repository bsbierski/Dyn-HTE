using JLD2,DelimitedFiles
include("../../plotConventions.jl")
include("../../LatticeGraphs.jl")
include("../../Embedding.jl")
include("../../ConvenienceFunctions.jl") 

### load graph evaluations and prepare lattice  
L = 12
n_max = 1*L
spin_length = 1/2
hte_graphs = load_dyn_hte_graphs(spin_length,L);

hte_lattice = getLattice(L,"square");
#display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

### compute all correlations in the lattice (or load them)
cd("CaseStudy/Square_Lattice_BSb/")
fileName_c = "Square_Lattice_c_iipDyn_nmax"*"_L"*string(L)*".jld2"
if isfile(fileName_c)
    c_iipDyn_mat = load_object(fileName_c)
else
    @time c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs);
    save_object(fileName_c,c_iipDyn_mat)
end


#########################################################################################
###### Dynamic structure factor (DSF) ###################################################
#########################################################################################

###### dynamical Matsubara correlator (k-space)
k,k_label = (1.0*π,0.0), "(pi,0)"
c_kDyn = get_c_kDyn(k,c_iipDyn_mat,hte_lattice)
m_vec = get_moments_from_c_kDyn(c_kDyn)
poly_x = Polynomial([0,1],:x)

x_vec_bare = collect(0:0.025:1.4)
x_vec = collect(0.0:0.1:4.0)

x0_vec = [1.0,2.0,3.0]      # for these x the DSF will be computed

### with u-series
if true
    r_max = 3                  # maximal order for moment
    
    f=0.7
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)
    u_vec = tanh.(f .* x_vec)
    u0_vec = tanh.(f .* x0_vec)
    m0_vec = [Float64[] for _ in x0_vec]

    plt_m = plot([0],[0],xlims=(0,x_vec[end]),ylims=(0,4.5),label="",xlabel=L"x=J/T",ylabel=L"x \cdot m_{\mathbf{k},2r}(x) \, / \, m_{\mathbf{k},2r}(0)",legend=:topleft);
    plot!(plt_m,x_vec,-x_vec,color=:grey,label="x bare");
    plot!(plt_m,x_vec,-x_vec,color=:grey,linestyle=linestyle_vec[2],label="u Padé [7-r,6-r]");
    plot!(plt_m,x_vec,-x_vec,color=:grey,linestyle=linestyle_vec[3],label="u Padé [6-r,5-r]");
    annotate!(plt_m,3,2,text(L"\mathbf{k}="*string(k_label)*",  f="*string(f),7));
    
    plot!(plt_m,title="SquareLattice AFM S=1/2: moment at k="*k_label*" (f=$f)")
    #plot!([2/f,2/f],[0,4],color=:grey,linewidth=10,label="u-freezing",alpha=0.3)
    for r in 0:r_max
        #xm_norm_r = m_vec[1+r]
        xm_norm_r = coeffs(poly_x * (m_vec[1+r]/m_vec[1+r](0)))
        println()
        println("r=$r")
        @show xm_norm_r
        p_u = Polynomial(ufromx_mat[1:n_max+2-2*r,1:n_max+2-2*r]*xm_norm_r)
        
        plot!(plt_m,x_vec_bare,Polynomial(xm_norm_r).(x_vec_bare),color=color_vec[r+1],linewidth=0.4,label="r=$r",alpha=0.7) 
        #plot!(plt_m,x_vec,p_u.(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[1],label="u-poly r=$r")
        
        plot!(plt_m,x_vec,get_pade(p_u,7-r,6-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="",alpha=0.7)
        plot!(plt_m,x_vec,get_pade(p_u,6-r,5-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="",alpha=0.7)

        ### extract moments at x0_vec 
        for x0_pos in eachindex(x0_vec)
            x0 = x0_vec[x0_pos]
            u0 = u0_vec[x0_pos]
            push!(m0_vec[x0_pos],m_vec[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
        end
    end

    xPlots,yPlots=1,1
    plt_final = plot(plt_m,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.65*aps_width*yPlots))
    display(plt_final)
    savefig(plt_final,"moments_u-series_k"*k_label*".png")
end


### δ_r and DSF for x ∈ x0_vec
if true
    w_vec = collect(0.0:0.025:5)
    plt_JS = plot([0,0],[0,0.01],color=:grey,label="Dyn-HTE",xlims=(0,w_vec[end]),xlabel=L"\omega/J=w",ylabel=L"J \, S(\mathbf{k}="*k_label*L",\omega)",legend=:topleft);#,title="SquareLat AFM S=1/2: "*"JS("*k_label*",w)")

    plt_δ =plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_{\mathbf{k},r}",legend=:topleft);

    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]

        ### plot Dyn-HTE
        δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
        scatter!(plt_δ,r_vec,δ_vec,color=thermalCol4_vec[x0_pos],label="x=$x0")
        δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,2000,true)
        plot!(plt_δ,r_max:7,δ_vec_ext[r_max+1:7+1],label="",color=thermalCol4_vec[x0_pos])
        plot!(plt_JS,w_vec,[JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec],color=thermalCol4_vec[x0_pos],label="")

        ### plot experimetnal data from DallaPiazza
        fileName = "DallaPiazza_exp_S(pi,0).csv"
        factor = 1200
        if isfile(fileName)
            data = readdlm(fileName,',',Float64)
            scatter!(plt_JS,data[:,1],data[:,2]/factor,color=:blue,marker=:cross,markersize=3.0,label="")#,label="x=$x0"
        end

    end
    xPlots,yPlots=1,2
    plt_final = plot(plt_δ, plt_JS,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
    display(plt_final)
    savefig(plt_final,"JS_k"*k_label*".png")
end


### DSF for k-path
if true
    w_vec = collect(0.0:0.025:3.5)
    r_max = 3                
    f=0.7
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)

    x0 = 2.0
    u0 = tanh.(f .* x0)

    ### define and generate k-path 
    path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
    pathticks = ["(π/2,π/2)","(π,0)","(π,π)","(π/2,π/2)","(0,0)","(π,0)"]

    Nk = 75  #75
    k_vec,kticks_positioins = create_brillouin_zone_path(path, Nk)
    JSkw_mat = zeros(Nk+1,length(w_vec))

    ### fill JSkw_mat
    for k_pos in eachindex(k_vec)
        k = k_vec[k_pos]

        c_kDyn = get_c_kDyn(k,c_iipDyn_mat,hte_lattice)
        m_vec = get_moments_from_c_kDyn(c_kDyn)
        m0 = Float64[]

        for r in 0:r_max
            xm_norm_r = coeffs(poly_x * (m_vec[1+r]/m_vec[1+r](0)))
            p_u = Polynomial(ufromx_mat[1:n_max+2-2*r,1:n_max+2-2*r]*xm_norm_r)
            
            push!(m0,m_vec[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
        end

        δ_vec,r_vec = fromMomentsToδ(m0)
        δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,2000,true)
        JSkw_mat[k_pos,:] = [JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec]

    end
end

### plot JS(k,ω)
if true
    using CairoMakie

    fig = Figure(fontsize=8,size=(aps_width,0.6*aps_width));
    ax=Axis(fig[1,1],xlabel=L"\mathbf{k}",ylabel=L"\omega/J=w",xlabelsize=8,ylabelsize=8);
    hm=CairoMakie.heatmap!(ax,collect(0:Nk)/(Nk),w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,0.4),highclip=:white);
    ax.xticks = ((kticks_positioins .- 1)/(Nk),pathticks)
    CairoMakie.Colorbar(fig[:, end+1], hm,size=8, label = L"J S(k,\omega)")
    CairoMakie.text!(ax,"x=J/T=$x0",position=[(0.05,0.5)],color=:white)
    CairoMakie.text!(ax,"f=$f",position=[(0.05,0.2)],color=:white)

    #Add experimental magnon data from DallaPiazza at T/J~ 0.09
    data = readdlm("DallaPiazza_exp_MagnonEnergy.csv",',',Float64)
    CairoMakie.plot!(ax,data[:,1],data[:,2],color=:purple,marker=:cross,markersize=5.0,label="")

    resize_to_layout!(fig);
    display(fig)

    save("Square_Lattice_JSkw.pdf",fig)
end