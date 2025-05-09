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

hte_lattice = getLattice(L,"triang");
#display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

### define special points in BZ
Γ,K,M = (0,0), (2*π/3,2*π/sqrt(3)), (0,2*π/sqrt(3))

### compute all correlations in the lattice (or load them)
fileName_c = "CaseStudy/Triangular_Lattice_BSb/Triangular_Lattice_c_iipDyn_nmax"*string(nmax)*"_L"*string(L)*".jld2"
if isfile(fileName_c)
    println("loading "*fileName_c)
    c_iipDyn_mat = load_object(fileName_c)
else
    @time c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs);
    save_object(fileName_c,c_iipDyn_mat)
end


#########################################################################################
###### Equal-time correlations  #########################################################
#########################################################################################
c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat)

if true #test uniform χ against HTE literature
    ###### HTE for χ = ∑_i<SiSj> from Dyn-HHTE ##############################################
    println( [sum(c_iipEqualTime_mat[j,1][n+1] for j in 1:hte_lattice.lattice.length) for n in 0:n_max]' )

    ### comparison to [Elstner,Singh,Young,PRL71.10(1993)] for S=1/2 only - ok [they have other definition for \chi and hence start at order β]
    coeffs_HTE = [1,-12,144,-1632,18000,-254016,5472096,-109168128, 818042112, 17982044160, 778741928448, -90462554542080, 829570427172864];
    println( [coeffs_HTE[1+n]*(-1)^n//4^(n+1)//factorial(Int128(n+1)) for n in 0:n_max]' )
end

###### equal-time correlations (k-space) and comparison to XTRG of [ChenPRB2019]
if true ### Gk for special k vs x with u-Pade
    k,k_label = K,"K"
    f=0.2  #f=0.3 for k=M, f=0.2 for k=K.
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    x_vec = collect(0:0.05:5.1)
    u_vec = tanh.(f .* x_vec)

    plt = plot(xlabel="x=J/T",legend=:bottomright,title="Triangular AFM S=1/2: equal-time struc-fac k=$k_label f=$f")

    coeffs_x = flipEvenIndexEntries(get_c_k(k , c_iipEqualTime_mat,hte_lattice))
    
    p_x = Polynomial(coeffs_x)
    p_u = Polynomial(ufromx_mat*coeffs_x)

    plot!(plt,x_vec[x_vec .< 1.15], Polynomial(coeffs_x).(x_vec[x_vec .< 1.15]),color=:blue,label="x-series",linewidth=0.5)
    plot!(plt,x_vec,get_pade(p_u,6,6).(u_vec),color=:blue,linestyle=linestyle_vec[2],label="u-Padé[6,6] (f=$f)")
    plot!(plt,x_vec,get_pade(p_u,5,5).(u_vec),color=:blue,linestyle=linestyle_vec[3],label="u-Padé[5,5] (f=$f)")
    #plot!(plt,x_vec,get_pade(p_u,4,4).(u_vec),color=:blue,linestyle=linestyle_vec[4],label="u-Padé[4,4]")
    #plot!(plt,x_vec,get_pade(p_u,6,5).(u_vec),color=:blue,linestyle=linestyle_vec[5],label="u-Padé[6,5]")

    plot!(plt,x_vec,get_pade(p_x,6,6).(x_vec),color=:green,linestyle=linestyle_vec[2],label="x-Padé[6,6]")
    plot!(plt,x_vec,get_pade(p_x,5,5).(x_vec),color=:green,linestyle=linestyle_vec[3],label="x-Padé[5,5]")
    #plot!(plt,x_vec,get_pade(p_x,4,4).(x_vec),color=:green,linestyle=linestyle_vec[4],label="x-Padé[4,4]")


    ### plot Gk from [ChenPRB2019]
    fileNameXTRG = "CaseStudy/Triangular_Lattice_BSb/Triangular_ThreeGk"*k_label*"_Chen.csv"
    if isfile(fileNameXTRG)
        data = readdlm(fileNameXTRG,',',Float64)
        scatter!(plt,1 ./ data[:,1],data[:,2] ./ 3,color=:blue,label="")
    end

    xPlots,yPlots=1,1
    plt_final = plot(plt,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
    display(plt_final)
    savefig(plt_final,"CaseStudy/Triangular_Lattice_BSb/Triangular_EqualTime_Gk"*k_label*".png")
end


#########################################################################################
###### Static structure factor (iν_m=0) #################################################
#########################################################################################
if true ### plot χT at k=K vs T
    k,k_label = K,"K"
    x_vec_bare = collect(0:0.01:1.2)
    x_vec = collect(0:0.01:3)


    plt_TχK_vs_x = plot(xlims=(0,x_vec[end]),ylims=(0.25,0.7),xlabel=L"x=J/T",ylabel=L"TG_{k=K}(i\nu_{m=0}) = Tχ_{k=K}",  label="",legend=:topleft)


    ### plot k=K data from [KulaginPRB2013]
    fileNameKulagin = "CaseStudy/Triangular_Lattice_BSb/Triangular_StaticSusc_Kulagin_atK.csv"
    if isfile(fileNameKulagin)
        data = readdlm(fileNameKulagin,',',Float64)
        scatter!(plt_TχK_vs_x, 1 ./ data[:,1], data[:,2] .* data[:,1] , color=:black,label="[Kulagin2013]")
    end

    ### x-Padé and u-Padé approximants
    m = 0
    p_x = sum([cos(dot(k,getSitePosition(hte_lattice.lattice,i).-getSitePosition(hte_lattice.lattice,hte_lattice.basis_positions[1])))*get_TGiip_Matsubara_xpoly(c_iipDyn_mat,i,1,m) for i in 1:hte_lattice.lattice.length])

    pade_vec = [[6,6],[7,5],[5,7],[5,5]]
    label=""
    for  (pade_pos,pade) in enumerate(pade_vec)
        label="x-Padé "*string(pade)
        y_vec = get_pade(p_x,pade[1],pade[2]).(x_vec)  
        plot!(plt_TχK_vs_x,x_vec,y_vec,color=:darkblue,alpha=0.6,linestyle=linestyle_vec[pade_pos],label=label)
    end
    for  (pade_pos,pade) in enumerate(pade_vec)
        f=0.25
        ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
        p_u = Polynomial(ufromx_mat*coeffs(p_x))
        u_vec = tanh.(f .* x_vec)
        if pade_pos==1 label="u-Padé (f=$f)" else label="" end
        plot!(plt_TχK_vs_x,x_vec,get_pade(p_u,pade[1],pade[2]).(u_vec),color=:purple,alpha=0.6,linestyle=linestyle_vec[pade_pos],label=label)
    end
    display(plt_TχK_vs_x)
end

if true
    ###### χ vs k #########
    x_vec = 1 ./ [2, 1, 0.5, 0.375]

    ### define and generate k-path 
    path = [Γ,K,M,Γ]
    pathticks = ["Γ","K","M","Γ"]

    Nk = 200  #75
    k_vec,kticks_positioins = create_brillouin_zone_path(path, Nk)

    plt_Tχk_vs_k=plot(xlims=(0,1),xlabel="",ylabel=L"Tχ_\mathbf{k} = TG_\mathbf{k}(i\nu_{m=0})", legendfontsize=5.4,  label="",legend=:topright,xticks=((kticks_positioins .- 1)/(Nk),pathticks))

    ### plot Dyn-HTE
    for (x_pos,x) in enumerate(x_vec)
        TGm0k_vec = zeros(length(k_vec))

        for (k_pos,k) in enumerate(k_vec)
            m = 0
            p_x = sum([cos(dot(k,getSitePosition(hte_lattice.lattice,i).-getSitePosition(hte_lattice.lattice,hte_lattice.basis_positions[1])))*get_TGiip_Matsubara_xpoly(c_iipDyn_mat,i,1,m) for i in 1:hte_lattice.lattice.length])

            pade_vec = [[6,6]]

            ### u-Padé approximants
            label=""
            for  (pade_pos,pade) in enumerate(pade_vec)
                
                f=0.25
                ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
                p_u = Polynomial(ufromx_mat*coeffs(p_x))
                u = tanh.(f .* x)

                TGm0k_vec[k_pos] = get_pade(p_u,pade[1],pade[2]).(u)
            end
        end

        plot!(plt_Tχk_vs_k,collect(0:Nk)/Nk,TGm0k_vec,color=thermalCol4_vec[x_pos],linealpha=0.75,label="T/J=1/x="*string(1/x))

        if true ### fit rMF
            function γtriang(kx,ky) 
                return 2 * (cos.(kx) .+ 2 .* cos.(kx ./ 2) .* cos.((3^(1/2) .* ky) ./ 2) )
            end
            
            kx_vec = [k[1] for k in k_vec]
            ky_vec = [k[2] for k in k_vec]

            function rMF(k_pos,p) 
                return 1 ./ (p[1] .+ p[2] .* γtriang(kx_vec[k_pos],ky_vec[k_pos])) 
            end 
            
            ### initial guess for params
            b1 = 1/2*(1/2+1)/3
            z = 6
            I3 =6
            f3rdOrder = 1/b1 + z*(1/6*(6*b1+1)*x^2+1/24*(4*b1+1)*x^3)- I3/6*b1*(6*b1+1)*x^3
            g3rdOrder = x + x^2 / 12 + x^3 /120 * (48*b1^2 + 16* b1 +3)
            
            p0 = [f3rdOrder, g3rdOrder]
            
            fit = LsqFit.curve_fit(rMF, collect(1:(Nk+1)), TGm0k_vec, p0)
            p=fit.param
            @show p

            plot!(collect(0:Nk)/Nk, rMF(collect(1:(Nk+1)),fit.param) ,color=thermalCol4_vec[x_pos], linestyle=:dash, label="rMF [f,g]="*string(round.(p,digits=2)) )
        end
    end

    if true ### plot bold-line diagMC data from [KulaginPRB2013]
        for (x_pos,x) in enumerate(x_vec)
            fileNameKulagin = "CaseStudy/Triangular_Lattice_BSb/Triangular_StaticSusc_Kulagin_GammaKMGamma_T"*string(1/x)*".csv"
            if isfile(fileNameKulagin)
                data = readdlm(fileNameKulagin,',',Float64)
                scatter!(plt_Tχk_vs_k, data[:,1], data[:,2] / x, markeralpha=0.4, color=thermalCol4_vec[x_pos],label="")
            end
        end
    end

    display(plt_Tχk_vs_k)
end


### prepare triangular lattice inset
if true
    plot!(plt_TχK_vs_x,inset=bbox(0.65,0.36,0.36,0.36),subplot=2)
    plot!(plt_TχK_vs_x[2],xlims=(-1.3,1.3),ylims=(-1.2,1.2),aspect_ratio = :equal,xaxis=false,yaxis=false)
    a1 = [1/2, sqrt(3)/2]
    a2 = [1, 0]
    a3 = [-1/2, sqrt(3)/2]
    for x in -5:1:5, y in -5:1:5
        r = x*a1 .+ y*a2
        r1 = r .+ a1
        r2 = r .+ a2
        r3 = r .+ a3
        plot!(plt_TχK_vs_x[2],[r[1],r1[1]],[r[2],r1[2]],markers=:dot,color=:grey,label="")
        plot!(plt_TχK_vs_x[2],[r[1],r2[1]],[r[2],r2[2]],markers=:dot,color=:grey,label="")
        plot!(plt_TχK_vs_x[2],[r[1],r3[1]],[r[2],r3[2]],markers=:dot,color=:grey,label="")
    end
end

### finalize
xPlots,yPlots=1,2
plt_final = plot(plt_TχK_vs_x,plt_Tχk_vs_k, layout=(yPlots,xPlots), size=(aps_width*xPlots,(0.63)*aps_width*yPlots),dpi=600)
display(plt_final)
savefig(plt_final,"CaseStudy/Triangular_Lattice_BSb/Triangular_StaticSF.png")




#########################################################################################
###### Dynamic structure factor (DSF) ###################################################
#########################################################################################

###### dynamical Matsubara correlator (k-space)
k,k_label = M,"M"
c_kDyn = get_c_k(k,c_iipDyn_mat,hte_lattice)
m_vec = get_moments_from_c_kDyn(c_kDyn)
poly_x = Polynomial([0,1],:x)

x_vec_bare = collect(0.0:0.025:1.2)
x_vec = collect(0.0:0.2:4.0)

x0_vec = 1 ./ [3.0,1.9,1.5,1.2,0.95,0.8,0.7,0.6,0.5,0.43,0.38]  # for these x the DSF will be computed

##### plot DSF and related quantities
if true
    w_vec = collect(0.0:0.02:3.7)

    ###### Pade for moments with x-series and u-series
    r_max = 3   # maximal order for moment
    f=0.55
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)
    u_vec = tanh.(f .* x_vec)
    u0_vec = tanh.(f .* x0_vec)
    m0_vec = [Float64[] for _ in x0_vec]

    plt_m = plot([0],[0],xlims=(0,x_vec[end]),label="")
    plot!(plt_m,xlabel=L"x=J/T",ylabel=L"x \cdot m_{\mathbf{k},2r}(x) \, / \, m_{\mathbf{k},2r}(0)",legend=:topleft)
    plot!(plt_m,-x_vec,0*x_vec,color=:grey,label="x bare")
    plot!(plt_m,-x_vec,0*x_vec,color=:grey,linestyle=linestyle_vec[2],label="u Padé [7-r,6-r]")
    plot!(plt_m,-x_vec,0*x_vec,color=:grey,linestyle=linestyle_vec[3],label="u Padé [6-r,5-r]")
    annotate!(plt_m,3,2,text(L"\mathbf{k}="*string(k_label)*",  f="*string(f),7))
    
    for r in 0:r_max
        println("r=$r")

        xm_norm_r = coeffs(poly_x * (m_vec[1+r]/m_vec[1+r](0)))
        p_x = Polynomial(xm_norm_r)
        p_u = Polynomial(ufromx_mat[1:n_max+2-2*r,1:n_max+2-2*r]*xm_norm_r)
        
        plot!(plt_m,x_vec_bare,Polynomial(xm_norm_r).(x_vec_bare),color=color_vec[r+1],linewidth=0.4,label="r=$r",alpha=0.7) 

        ### x-Padé moments are not well behaved (grow large or negative)
        #plot!(plt_m,x_vec,get_pade(p_x,7-r,6-r).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="",alpha=0.7)
        #plot!(plt_m,x_vec,get_pade(p_x,6-r,5-r).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="",alpha=0.7)

        plot!(plt_m,x_vec,get_pade(p_u,7-r,6-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="",alpha=0.7)
        plot!(plt_m,x_vec,get_pade(p_u,6-r,5-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="",alpha=0.7)

        ### extract moments from u-Padé [7-r,6-r] at x0_vec 
        for x0_pos in eachindex(x0_vec)
            x0 = x0_vec[x0_pos]
            u0 = u0_vec[x0_pos]
            push!(m0_vec[x0_pos],m_vec[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
        end
    end

    ###### δ_r, JS and A for x ∈ x0_vec
    plt_δ=plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_{\mathbf{k},r}",legend=:bottomright)
    plt_JS = plot(xlims=(0,w_vec[end]),xlabel=L"\omega/J=w",ylabel=L"J \, S(\mathbf{k}="*k_label*L",\omega)",legendfontsize=5.0,legend=:topright)
    plt_JAo2π = plot(xlabel=L"\omega/J=w",ylabel=L"J \, A(\mathbf{k}="*k_label*L",\omega)",legend=:topleft)

    #plt_JSw0=plot([0.55,0.55],[0.0,0.09],xscale=:log10,xlims=(0.1,1.02/x0_vec[1]),label="roton-like energy [Zheng2006]",color=:grey,xlabel=L"T/J=1/x",ylabel=L"J \, S(\mathbf{k}="*k_label*L",\omega \rightarrow 0)",legend=:bottomright)

    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]

        ### plot Dyn-HTE
        δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
        scatter!(plt_δ,r_vec,δ_vec,color=thermalCol13_vec[x0_pos],label="")
        δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,4000,true)
        plot!(plt_δ,r_max+1:6,δ_vec_ext[r_max+2:7],label="",color=thermalCol13_vec[x0_pos])

        JSw_vec = [JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec]

        if k_label=="K"
            plot!(plt_JS,w_vec, JSw_vec,color=thermalCol13_vec[x0_pos],label="")
        else
            plot!(plt_JS,w_vec, JSw_vec,color=thermalCol13_vec[x0_pos],label="T/J="*string(round(1/x0,digits=2)))
        end

        plot!(plt_JAo2π,w_vec, JSw_vec .* (1 .- exp.(-x0 .* w_vec)),color=thermalCol13_vec[x0_pos],label="")
        #scatter!(plt_JSw0,[1/x0],[JSw_vec[1]],color=thermalCol15_vec[x0_pos],label="")
    end

    ### plot results
    xPlots,yPlots=4,1
    plt_final = plot(plt_m, plt_δ, plt_JS, plt_JAo2π, layout=(yPlots,xPlots), size=(aps_width*xPlots,0.7*aps_width*yPlots))
    display(plt_final)
    savefig(plt_final,"CaseStudy/Triangular_Lattice_BSb/Triangular_DSF_k"*k_label*".png")
end

### prepare insets
if k_label=="K" 
    plot!(plt_JS,inset=bbox(0.62,0.03,0.36,0.49),subplot=2) 
end
if k_label=="M" 
    plot!(plt_JS,xlabel="",xformatter=_->"") 
    plot!(plt_JS,inset=bbox(0.2,0.65,0.3,0.3),subplot=2)
    plot!(plt_JS[2],xlims=(-1.3,1.3),ylims=(-1.2,1.2),aspect_ratio = :equal,xaxis=false,yaxis=false)
    a1 = [1/2, sqrt(3)/2]
    a2 = [1, 0]
    a3 = [-1/2, sqrt(3)/2]
    for x in -5:1:5, y in -5:1:5
        r = x*a1 .+ y*a2
        r1 = r .+ a1
        r2 = r .+ a2
        r3 = r .+ a3
        plot!(plt_JS[2],[r[1],r1[1]],[r[2],r1[2]],markers=:dot,color=:black,label="")
        plot!(plt_JS[2],[r[1],r2[1]],[r[2],r2[2]],markers=:dot,color=:black,label="")
        plot!(plt_JS[2],[r[1],r3[1]],[r[2],r3[2]],markers=:dot,color=:black,label="")
    end
end

### scaling plot of DSF at k=K (as inset)
if k_label=="K" && true
    w_max = 1.0
    w_vec = collect(0.0:0.02:w_max)
    α = 1.1

    annotate!(plt_JS[2],1,0.1,text(L"\alpha="*string(α)*"0(2)",7))
    plot!(plt_JS[2],xlims=(0,2.7),ylims=(0.08,0.165),xlabel=L"\omega/T\;\;(\omega \leq"*string(w_max)*L" \,J)",ylabel=L"J\, S(\mathbf{k}="*k_label*L",\omega) \cdot (T/J)^\alpha")

    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]
        if 1/x0 >= 0.3 && 1/x0 <= 1.0

            ### plot Dyn-HTE
            δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
            scatter!(plt_δ,r_vec,δ_vec,color=thermalCol13_vec[x0_pos],label="x=$x0")
            δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,4000,true)
            JSw_vec = [JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec]

            plot!(plt_JS[2],w_vec .* x0, JSw_vec .* (1/x0)^α ,color=thermalCol13_vec[x0_pos],label="")
        end
    end
    display(plt_JS)
end

### final plot for paper
if k_label=="K" plt_JS_K = deepcopy(plt_JS) end
if k_label=="M" plt_JS_M = deepcopy(plt_JS) end


###### run the above for K and M and then put together
xPlots,yPlots=1,2
plt_final = plot(plt_JS_M,plt_JS_K, layout=(yPlots,xPlots), size=(aps_width*xPlots,(0.48)*aps_width*yPlots),dpi=600)
display(plt_final)
savefig(plt_final,"CaseStudy/Triangular_Lattice_BSb/Triangular_DSF.png")




### background information: scaling plot of DSF at k=K at three different α (in quantum critical T-regime?)
if k_label=="K" && false
    w_max = 1.0
    w_vec = collect(0.0:0.02:w_max)
    α1,α2,α3 = [1.0,1.1,1.2]

    plt_JS_scaled1 = plot(ylims=(0.08,0.17),legend=:topright)
    plt_JS_scaled2 = plot(ylims=(0.08,0.17),legend=:topright)
    plt_JS_scaled3 = plot(ylims=(0.08,0.17),legend=:topright)
    plot!(plt_JS_scaled1,xlabel=L"\omega/T\;\;(\omega \leq"*string(w_max)*L" \,J)",ylabel=L"J\, S(\mathbf{k}="*k_label*L",\omega) \times (T/J)^\alpha  \;\;\; \alpha="*string(α1))
    plot!(plt_JS_scaled2,xlabel=L"\omega/T\;\;(\omega \leq"*string(w_max)*L" \,J)",ylabel=L"J\, S(\mathbf{k}="*k_label*L",\omega) \times (T/J)^\alpha  \;\;\; \alpha="*string(α2))
    plot!(plt_JS_scaled3,xlabel=L"\omega/T\;\;(\omega \leq"*string(w_max)*L" \,J)",ylabel=L"J\, S(\mathbf{k}="*k_label*L",\omega) \times (T/J)^\alpha  \;\;\; \alpha="*string(α3))


    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]
        if 1/x0 >= 0.3 && 1/x0 <= 1.0

            ### plot Dyn-HTE
            δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
            scatter!(plt_δ,r_vec,δ_vec,color=thermalCol13_vec[x0_pos],label="x=$x0")
            δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,4000,true)
            JSw_vec = [JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec]

            plot!(plt_JS_scaled1,w_vec .* x0, JSw_vec .* (1/x0)^α1 ,color=thermalCol13_vec[x0_pos],label="T/J="*string(round(1/x0,digits=2)))
            plot!(plt_JS_scaled2,w_vec .* x0, JSw_vec .* (1/x0)^α2 ,color=thermalCol13_vec[x0_pos],label="T/J="*string(round(1/x0,digits=2)))
            plot!(plt_JS_scaled3,w_vec .* x0, JSw_vec .* (1/x0)^α3 ,color=thermalCol13_vec[x0_pos],label="T/J="*string(round(1/x0,digits=2)))

        end
    end

    xPlots,yPlots=3,1
    plt_final = plot(plt_JS_scaled1,plt_JS_scaled2,plt_JS_scaled3, layout=(yPlots,xPlots), size=(aps_width*xPlots,0.5*aps_width*yPlots))
    display(plt_final)
    savefig(plt_final,"CaseStudy/Triangular_Lattice_BSb/Triangular_JS_k"*k_label*"_scaling.png")
end


