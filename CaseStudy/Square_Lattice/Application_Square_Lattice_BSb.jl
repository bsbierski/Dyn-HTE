using JLD2,DelimitedFiles
include("../../plotConventions.jl")
include("../../Embedding.jl")
include("../../LatticeGraphs.jl")
include("../../ConvenienceFunctions.jl") 

### maximum expansion order
n_max = 12                                      

### load list of unique graphs
gG_vec_unique = give_unique_gG_vec(n_max);

### load dictionaries of all lower orders C_Dict_vec 
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,n_max+1);
for ord = 0:n_max
    C_Dict_vec[ord+1]  = load_object("GraphEvaluations/Spin_S1half/C_"*string(ord)*".jld2")
end 

### define lattice ball for embedding (it is enough for embedding of n_max graphs to have ball radius L=n_max)
L = 12
lattice,LatGraph,center_sites = getLattice_Ball(L,"square");
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.11,fontsize=4,nodeshape=:rect,curves=false))

### compute all correlations in the lattice (or load them)
fileName_c = "CaseStudy/Square_Lattice/Square_c_iipDyn_nmax"*string(n_max)*"_L"*string(L)*".jld2"
if isfile(fileName_c)
    c_iipDyn_mat = load_object(fileName_c)
else
    c_iipDyn_mat = get_c_iipDyn_mat(LatGraph,lattice,center_sites,n_max,gG_vec_unique,C_Dict_vec);
    save_object("CaseStudy/Square_Lattice/Square_c_iipDyn.jld2",c_iipDyn_mat)
end

### define special points in BZ
Γ,X,M = [0,0], [π,0], [π,π]
k_vec,klabel_vec = [Γ,X,M],["Γ","X","M"]

#########################################################################################
###### Equal-time correlations  #########################################################
#########################################################################################
### TODO: Adapt all of this
c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat,n_max)

if false #test uniform χ against HTE literature
    ###### HTE for χ = ∑_i<SiSj> from Dyn-HHTE ##############################################
    println( [sum(c_iipEqualTime_mat[:,1,n+1]) for n in 0:n_max]' )

    # comparison to [Lohmann,PRB,89.014415(2014)] for S=1/2 for n<=9 - ok
    function cLohmann(nn)::Rational{Int}
        n = nn+1 #Lohmann uses a different convention for order
        r= 3//4
        if n == 1
            return (1//3) * r
        elseif n == 2
            return -(4//9) * r^2
        elseif n == 3
            return (1//9) * r^2 * (-1 + 4r)
        elseif n == 4
            return -(4//405) * r^2 * (3 - 28r + 37r^2)
        elseif n == 5
            return (1//4860) * r^2 * (-45 + 702r - 1892r^2 + 1328r^3)
        elseif n == 6
            return -(1//510300) * r^2 * (1728 - 35946r + 164289r^2 - 207896r^3 + 102576r^4)
        elseif n == 7
            return (1//6123600) * r^2 * (-8694 + 218916r - 1401381r^2 + 2888772r^3 - 2251248r^4 + 909184r^5)
        elseif n == 8
            return -(1//22963500) * r^2 * (15390 - 446256r + 3538764r^2 - 10535337r^3 + 12202552r^4 - 7318640r^5 + 2416640r^6)
        elseif n == 9
            return (1//7715736000) * r^2 * (-2710665 + 87954822r - 807482331r^2 + 3091042674r^3 - 5118502560r^4 + 4009481184r^5 - 2113197952r^6 + 518354176r^7)
        elseif n == 10
            return −1//1273096440000*r^2 * (257596200 − 9180862110r + 93799827171r^2 − 426255134022r^3 + 931126345494r^4 − 977085756168r^5 + 621427831616r^6 − 280517703040r^7 + 48779713280r^8)
        else
            return 0  # Undefined for n > 10 in given series
        end
    end
    println(  [cLohmann(n)*(-1)^n for n in 0:9]' )

    ## comparison to [Elstner and Young,PRB50.6871(1994)] for S=1/2 only - ok
    ## [they have other definition for \chi and hence start at order β. They have a typo (fixed here)]
    coeffs = [4,-32,192,-384,-1280,-155136,2711296,56705024,-1716811776,-47711784960,2004747075584,55843726884864,-3367208347123712];
    println( [coeffs[1+n]*(-1)^n//4^(n+2)//factorial(n+1) for n in 0:n_max]' )
end

##### equal-time correlations (r-space)
if false #standard Padé in x
    j_vec = [191,154,153,151]
    x_vec_bare = collect(0:0.025:1.35)
    x_vec = collect(0:0.1:50)
    plt = plot([0],[0],label="",xlabel="x=J/T",ylims=(-0.07,0.05),legend=:right,title="Kagome AFM S=1/2: equal-time struc-fac L$L")
    for (j_pos,j) in enumerate(j_vec)

        @show j
        GjEqualTime_poly_x = Polynomial(1.0*flipEvenIndexEntries(c_iipEqualTime_mat[j,1,:]))
        @show flipEvenIndexEntries(c_iipEqualTime_mat[j,1,:])

        #plot!(plt,x_vec_bare,GjEqualTime_poly_x.(x_vec_bare),color=color_vec[j_pos],label="j=$j",lw=0.4)
        plot!(plt,x_vec,get_pade(GjEqualTime_poly_x,4,4).(x_vec),color=color_vec[j_pos],linestyle=linestyle_vec[2],label="[4,4] j=$j")
        plot!(plt,x_vec,get_pade(GjEqualTime_poly_x,6,6).(x_vec),color=color_vec[j_pos],linestyle=linestyle_vec[4],label="[6,6] j=$j")

        #plot!(plt,x_vec,get_intDiffApprox(GjEqualTime_poly_x, x_vec, 3,3,3),color=color_vec[j_pos],linestyle=linestyle_vec[5],label="")
    end
    xPlots,yPlots=1,1
    plt_final = plot(plt,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
    display(plt_final)
    #savefig(plt_final,"CaseStudy/Kagome/KagomeEqualTimeCorrelator_j_xsweep_L$L.png")
end

if false  # u=tanh(fx) with bare series in u or Pade in u
    j_vec = [191,153,154,151,155,120]

    f=0.72
    @variables x u
    x = taylor(atanh(u)/f, u, 0:12, rationalize=false)

    u_vec = 1 .- 10 .^ range(-6,stop=0,length=50)

    plt = plot([0],[0],label="",xlabel="x=J/T (f=$f)",ylims=(-0.07,0.03),xlims=(0,10.0),legend=:right,title="Kagome AFM S=1/2: equal-time corr L$L, Padé in u=tanh(fx)")
    for (j_pos,j) in enumerate(j_vec)

        @show j
        coeffs_x = 1.0*flipEvenIndexEntries(c_iipEqualTime_mat[j,1,:])

        p_u_ext = simplify(series(coeffs_x,x);expand=true)
        p_u = Polynomial(Symbolics.value.(taylor_coeff(p_u_ext,u,0:12,rationalize=false)),:u)

        plot!(plt,atanh.(u_vec)/f, p_u.(u_vec),color=color_vec[j_pos],linestyle=linestyle_vec[1],label="j=$j")

        plot!(plt,atanh.(u_vec)/f, get_pade(p_u,4,4).(u_vec),color=color_vec[j_pos],linestyle=linestyle_vec[2],label="")
        plot!(plt,atanh.(u_vec)/f, get_pade(p_u,5,7).(u_vec),color=color_vec[j_pos],linestyle=linestyle_vec[3],label="")
        plot!(plt,atanh.(u_vec)/f, get_pade(p_u,6,6).(u_vec),color=color_vec[j_pos],linestyle=linestyle_vec[4],label="")

    end
    xPlots,yPlots=1,1
    plt_final = plot(plt,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
    display(plt_final)
    savefig(plt_final,"CaseStudy/Kagome/KagomeEqualTimeCorrelator_j_usweep_f"*string(f)*"_L$L.png")
end

###### equal-time correlations (k-space)
if false ### Gk for special k vs x with u-Pade
    f=0.25
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    x_vec = collect(0:0.1:10)
    u_vec = tanh.(f .* x_vec)

    plt = plot([0],[0],ylims=(0,0.5),label="",xlabel="x=J/T",legend=:right,title="Kagome AFM S=1/2: equal-time struc-fac f=$f")
    for (k_pos,k) in enumerate(k_vec)

        coeffs_x = flipEvenIndexEntries(get_c_kEqualTime(k , c_iipEqualTime_mat,lattice,center_sites))
        p_u = Polynomial(ufromx_mat*coeffs_x)

        plot!(plt,x_vec, p_u.(u_vec),color=color_vec[k_pos],label="k="*klabel_vec[k_pos])
        plot!(plt,x_vec,get_pade(p_u,4,4).(u_vec),color=color_vec[k_pos],linestyle=linestyle_vec[2],label="")
        plot!(plt,x_vec,get_pade(p_u,6,6).(u_vec),color=color_vec[k_pos],linestyle=linestyle_vec[4],label="")

        #plot!(plt,x_vec,get_intDiffApprox(GkEqualTime_poly_x, x_vec, 3,3,3),color=color_vec[k_pos],linestyle=linestyle_vec[5],label="")
    end
    xPlots,yPlots=1,1
    plt_final = plot(plt,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
    display(plt_final)
    #savefig(plt_final,"CaseStudy/Kagome/Kagome_EqualTime_Gk_f$f"*".png")
end

if false ### BZ Gk plot
    f=0.47
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    x = 4.0
    u = tanh.(f .* x)
    pade_order=[6,6]

    LL=60
    step_vec = collect(-3*LL:1:3*LL)/LL
    function kk(f1,f2) 
        #f1 .* G1 + f2 .* G2
        f1 .* [π,0] + f2 .* [0,π] 
    end
    k_vec = [kk(f1,f2) for f1 in step_vec for f2 in step_vec if abs(kk(f1,f2)[1]/π)<=1.5 && abs(kk(f1,f2)[2]/π)<=1.5]
    kx_vec, ky_vec = [k[1] for k in k_vec],[k[2] for k in k_vec]

    ### equal-time correlator
    coeffs_x_vec = [flipEvenIndexEntries(get_c_kEqualTime(k,c_iipEqualTime_mat,lattice,center_sites)) for k in k_vec]

    if f==0.0 ### use x-Polynomial
        z_data = [get_pade(Polynomial(coeffs_x),pade_order...)(x) for coeffs_x in coeffs_x_vec]
    else ### use u-Polynomial
        z_data = [get_pade(Polynomial(ufromx_mat*coeffs_x),pade_order...)(u) for coeffs_x in coeffs_x_vec]
        #z_data = [Polynomial(ufromx_mat*coeffs_x)(u) for coeffs_x in coeffs_x_vec]
    end

    ### compute sum-rule
    Σ = 0.0
    VBZ = 6 * 4π/3 * π/3^0.5
    for (k_pos,k) in enumerate(k_vec)
        if k[1] >=0 && k[1] <= 4/3*π
            if k[2] <= 2/3^0.5 / (2/3)*k[1] && k[2]>=0.0 && k[2] <= 4*π/3^0.5 - k[1] * 4*π/3^0.5 / (4*π/3) 
                Σ += 6/VBZ * z_data[k_pos] * (π/LL)^2
            end
        end
    end

    ### prepare plot
    using CairoMakie
    xPlots,yPlots=1,1
    fig = Figure(size=(1*aps_width*xPlots,1*0.8*aps_width*yPlots),fontsize=8)
    ax=Axis(fig[1,1],xlabel=L"k_x/\pi",ylabel=L"k_y/\pi",xlabelsize=9,ylabelsize=9,aspect=1,title="Kagome: Padé$pade_order x="*string(x)*" f=$f"*" Sigma="*string(round(Σ,digits=3)))
    hm=CairoMakie.heatmap!(ax,kx_vec/π,ky_vec/π,z_data,colormap=:viridis,colorrange=(0.0,0.44),highclip=:white,lowclip=:black)
    lines!(ax,[4/3*cos(α*π/3) for α in 0:6],[4/3*sin(α*π/3) for α in 0:6],color=:grey)
    CairoMakie.Colorbar(fig[:, end+1], hm, label=L"\mathrm{equal-time \;\; structure \;\; factor}\;\;G_\mathbf{k}",)
    text!(ax,M[1]/π-0.1,M[2]/π-0.2;text=L"M",fontsize=12)
    text!(ax,K[1]/π-0.22,K[2]/π-0.1;text=L"K",fontsize=12)
    CairoMakie.ylims!(ax,(-1.5,1.5))
    CairoMakie.xlims!(ax,(-1.5,1.5))
    resize_to_layout!(fig)
    display(fig)
    save("CaseStudy/Kagome/Kagome_EqualTimeGk_pade"*string(pade_order)*"_f$f"*"_x$x.png",fig)
end

#########################################################################################
###### Dynamic structure factor (DSF) ###################################################
#########################################################################################

###### dynamical Matsubara correlator (k-space)
k,k_label = X,"X"
c_kDyn = get_c_kDyn(k,c_iipDyn_mat,lattice,center_sites)
m_vec = get_moments_from_c_kDyn_mat(c_kDyn)
poly_x = Polynomial([0,1],:x)

x_vec_bare = collect(0:0.025:1.4)
x_vec = collect(0.0:0.1:4.0)

### with x-series
if true
    plt_m = plot([0],[0],xlims=(0,x_vec[end]),ylims=(0,4),label="",xlabel="x=J/T",ylabel=L"x \cdot m_r(x)/m_r(0)",legend=:topleft,title="SquareLat AFM S=1/2: moments @k="*k_label)
    for r in 0:3
        #xm_norm_r = m_vec[1+r]
        xm_norm_r = (m_vec[1+r]/m_vec[1+r](0)) * poly_x
        println()
        println("r=$r")
        @show xm_norm_r
        plot!(plt_m,x_vec_bare,xm_norm_r.(x_vec_bare),color=color_vec[r+1],linewidth=0.4,label="r=$r")
        
        plot!(plt_m,x_vec,get_pade(xm_norm_r,7-r,6-r).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="")
        plot!(plt_m,x_vec,get_pade(xm_norm_r,6-r,5-r).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="")
    end
    xPlots,yPlots=1,1
    plt_final = plot(plt_m,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
    display(plt_final)
    #savefig(plt_final,"CaseStudy/Square_Lattice/moments_x-series_k"*k_label*".png")
end

### with u-series
if true
    x0_vec = [1.0,1.5,2.0,2.5]  # for these x the DSF will be computed
    r_max = 2                   # maximal order for moment
    
    f=0.4
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)
    u_vec = tanh.(f .* x_vec)
    u0_vec = tanh.(f .* x0_vec)
    m0_vec = [Float64[] for _ in x0_vec]

    plt_m = plot([0],[0],xlims=(0,x_vec[end]),ylims=(0,5),label="",xlabel=L"x=J/T",ylabel=L"x \cdot m_{\mathbf{k},2r}(x) \, / \, m_{\mathbf{k},2r}(0)",legend=:topleft);
    plot!(plt_m,x_vec,-x_vec,color=:grey,label="x bare");
    plot!(plt_m,x_vec,-x_vec,color=:grey,linestyle=linestyle_vec[2],label="u Padé [7-r,6-r]");
    plot!(plt_m,x_vec,-x_vec,color=:grey,linestyle=linestyle_vec[3],label="u Padé [6-r,5-r]");
    plot!(plt_m,x_vec,-x_vec,color=:grey,linestyle=linestyle_vec[3],label="u Padé [5-r,4-r]");
    annotate!(plt_m,3,2,text(L"\mathbf{k}="*string(k_label)*",  f="*string(f),7));
    
    #plot!(plt_m,title="Kagome AFM S=1/2: moment at k="*k_label*" (f=$f)")
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
        plot!(plt_m,x_vec,get_pade(p_u,5-r,4-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[4],label="",alpha=0.7)

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
    #savefig(plt_final,"CaseStudy/Kagome/Kagome_moments_u-series_k"*k_label*".png")
end

### DSF and δ_r (as inset) for x ∈ x0_vec
if true
    w_vec = collect(0.0:0.05:5)
    plt_JS = plot([0,0],[0,0.01],color=:grey,label="Dyn-HTE",xlims=(0,w_vec[end]),xlabel=L"\omega/J=w",ylabel=L"J \, S(\mathbf{k}=M,\omega)",legend=:topleft);#,title="SquareLat AFM S=1/2: "*"JS("*k_label*",w)")
    plt_δ =plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_{\mathbf{k},r}",legend=:topleft);
    #scatter!(plt_JS,[0],[-2],color=:grey,marker=:cross,label="NLCE+Gauss [Sherman2018]")

    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]

        ### plot Dyn-HTE
        δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
        scatter!(plt_δ,r_vec,δ_vec,color=thermalCol4_vec[x0_pos],label="x=$x0")
        δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,2000,true)
        plot!(plt_δ,r_max:7,δ_vec_ext[r_max+1:7+1],label="",color=thermalCol4_vec[x0_pos])
        plot!(plt_JS,w_vec,[JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec],color=thermalCol4_vec[x0_pos],label="")


        ### plot Sherman's NLCE
        #fileNameNCLE = "CaseStudy/Kagome/Sherman2018NLCE_data/Sherman_NLCE_JS_x"*string(x0)*"_"*k_label*".csv"
        #if isfile(fileNameNCLE)
        #    Sherman = readdlm(fileNameNCLE,',',Float64)
        #    scatter!(plt_JS,Sherman[:,1],Sherman[:,2],color=thermalCol4_vec[x0_pos],marker=:cross,markersize=3.0,label="")#,label="x=$x0"
        #end

    end
    xPlots,yPlots=1,2
    plt_final = plot(plt_δ, plt_JS,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
    display(plt_final)
    #savefig(plt_final,"CaseStudy/Kagome/Kagome_JS_k"*k_label*".png")
end

### plot for paper
#addABC(plt_m,"(a)")
#addABC(plt_JS[1],"(b)")
#xPlots,yPlots=1,2
#plt_final = plot(plt_m,plt_JS,  layout=Plots.grid(2,1,heights=[0.45,0.55]), size=(aps_width*xPlots,(0.45+0.62)/2*aps_width*yPlots))
#display(plt_final)
#savefig(plt_final,"CaseStudy/Kagome/Kagome_JS_k"*k_label*"_paper.pdf")