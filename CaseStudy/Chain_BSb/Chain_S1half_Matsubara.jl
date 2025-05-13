using JLD2,DelimitedFiles, HDF5
include("../../plotConventions.jl")
include("../../LatticeGraphs.jl")
include("../../Embedding.jl")
include("../../ConvenienceFunctions.jl") 

### load graph evaluations and prepare lattice  
L = 12
n_max = 1*L
spin_length = 1/2
hte_graphs = load_dyn_hte_graphs(spin_length,L);

hte_lattice = getLattice(L,"chain");
#display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

### compute all correlations in the lattice
c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs)

#########################################################################################
###### Matsubara correlator at m=0,1,2 in real-space and vs QMC
#########################################################################################

###### load QMC data at m=0,1,2 (have x=βJ=[1,2,...,8]; but one can use less) 
x_vec_QMC = collect(1:8)
imax=16
m_vec = [0,1,2]
TGiip_m_QMC = zeros((length(x_vec_QMC),imax+1,length(m_vec)))
TGiip_m_QMC_err = zeros((length(x_vec_QMC),imax+1,length(m_vec)))

for x_QMC_pos in eachindex(x_vec_QMC) 
    x_QMC=x_vec_QMC[x_QMC_pos]
    fid = h5open("CaseStudy/Chain_BSb/QMC_Worm/SpinHalfAFMHeisenbergChain/job_BSb_beta"*string(x_QMC)*".out.h5", "r")
    for (m_pos,m) in enumerate(m_vec)
        TGiip_m_QMC[x_QMC_pos,:,m_pos] = read(fid["simulation"]["results"]["DensDens_CorrFun_w$m"]["mean"]["value"])[1:imax+1]
        TGiip_m_QMC_err[x_QMC_pos,:,m_pos] = read(fid["simulation"]["results"]["DensDens_CorrFun_w$m"]["mean"]["error"])[1:imax+1]
    end
end

###### plot QMC data against Dyn-HTE
x_vec = collect(0:0.05:8.15)

###### plot Matsubara correlator in r-space
if true
    ## real space i-i'=di
    di_vec = [0,1,2]
    x_vec_bare = collect(0:0.01:1.93)
    
    ### prepare panels for m=0,1,2
    plt_m0=plot(x_vec,0*x_vec,color=:black,lw=0.8,xlims=(0,x_vec[end]),ylims=(-0.11,0.26), xformatter=:none,ylabel=L"TG_{ii^\prime}(i\nu_{m=0})",  label="",legend=:topright)
    plt_m1=plot(x_vec,0*x_vec,color=:black,lw=0.8,xlims=(0,x_vec[end]),ylims=(-0.035,0.04),xformatter=:none,ylabel=L"TG_{ii^\prime}({i\nu_{m=1}})",label="",legend=:bottomleft)
    plt_m2=plot(x_vec,0*x_vec,color=:black,lw=0.8,xlims=(0,x_vec[end]),ylims=(-0.01,0.02),xlabel=L"x=J/T",ylabel=L"TG_{ii^\prime}(i\nu_{m=2})",   label="",legend=:topleft)
    plt_vec = [plt_m0,plt_m1,plt_m2]

    ### plot Dyn-HTE results
    for di in di_vec, (m_pos,m) in enumerate(m_vec)
        
        p_x = get_TGiip_Matsubara_xpoly(c_iipDyn_mat,hte_lattice.basis_positions[1]+di,1,m)

        ### plot bare truncated series
        label=""
        if (di==0 && m==0)  label = L"\mathrm{bare \; series} \; n_{max}="*string(n_max) end    
        plot!(plt_vec[m_pos],x_vec_bare,p_x.(x_vec_bare),label=label,color=:grey,linewidth=0.5)

        pade_vec = [[6,6],[5,5]]

        ### x-Padé approximants
        label=""
        for  (pade_pos,pade) in enumerate(pade_vec)
            if di==0 && m==1 label="x-Padé "*string(pade) end
            y_vec = get_pade(p_x,pade[1],pade[2]).(x_vec)

            x_vec_show = 1*x_vec
            y_vec_show = 1*y_vec
            for p in 1:length(y_vec)
                if abs(y_vec_show[p])>0.251
                    x_vec_show = 1*x_vec_show[1:p]
                    y_vec_show = 1*y_vec_show[1:p]
                    break
                end
            end
            
            plot!(plt_vec[m_pos],x_vec_show,y_vec_show,color=:darkblue,alpha=0.6,linestyle=linestyle_vec[pade_pos],label=label)
        end

        ### u-Padé approximants
        f=0.205
        ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
        p_u = Polynomial(ufromx_mat*coeffs(p_x))
        u_vec = tanh.(f .* x_vec)

        label=""
        for  (pade_pos,pade) in enumerate(pade_vec)
            if di==0 && m==2 label="u-Padé "*string(pade)*" (f=$f)" end
            plot!(plt_vec[m_pos],x_vec,get_pade(p_u,pade[1],pade[2]).(u_vec),color=:purple,alpha=0.6,linestyle=linestyle_vec[pade_pos],label=label)
        end
    end

    ### plot QMC data
    for di in di_vec, (m_pos,m) in enumerate(m_vec)
        scatter!(plt_vec[m_pos],x_vec_QMC,TGiip_m_QMC[:,1+di,1+m],color=color_vec[di+1],label="",markersize=6)
        if m==0 annotate!(plt_vec[m_pos],x_vec_QMC[3]+0.5,TGiip_m_QMC[3,1+di,1+m]+0.02,text("i-i'="*string(di),7,color_vec[di+1])) end
    end

    ### integrated-differential approximants
    if false
        for di in di_vec, (m_pos,m) in enumerate(m_vec[1])

            p = get_TGiip_m_bare(c_iipDyn_mat,m,12)[center_sites[1]+di,1]
            
            for (IDA_pos,IDA) in enumerate([[2,4,4]])
                label=""
                if (m_pos==1 && di==0) label="IDA "*string(IDA) end
                plot!(plt_vec[m_pos],x_vec,get_intDiffApprox(p,x_vec,IDA...),linestyle=linestyle_vec[IDA_pos],label=label,color=:magenta,alpha=0.5)
            end
        end
    end

    ### put panels together
    xPlots,yPlots=1,3
    plt_final = plot(plt_vec...,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.4*aps_width*yPlots), dpi=600)
    display(plt_final)
    savefig(plt_final,"CaseStudy/Chain_BSb/Chain_S1half_Matsubara_r-space.png")
end

###### plot Matsubara correlator in k-space
if true
    ## momenta
    koπ_vec = [0.4, 1.0]
    x_vec_bare = collect(0:0.01:1.96)
    
    ### prepare panels for m=0,1,2
    plt_m0=plot(xlims=(0,x_vec[end]),ylims=(0,0.75),xformatter=:none,ylabel=L"TG_{k}(i\nu_{m=0})",  label="",legend=:topleft)
    plt_m1=plot(xlims=(0,x_vec[end]),ylims=(0,0.10),xformatter=:none,ylabel=L"TG_{k}({i\nu_{m=1}})",label="",legend=:topleft)
    plt_m2=plot(xlims=(0,x_vec[end]),ylims=(0,0.04),xlabel=L"x=J/T",ylabel=L"TG_{k}(i\nu_{m=2})",   label="",legend=:topleft)
    plt_vec = [plt_m0,plt_m1,plt_m2]

    ### plot Dyn-HTE results
    for (koπ_pos,koπ) in enumerate(koπ_vec), (m_pos,m) in enumerate(m_vec)
        
        b=hte_lattice.basis_positions[1]
        p_x = sum([cos(di*koπ*π)*get_TGiip_Matsubara_xpoly(c_iipDyn_mat,b+di,1,m) for di in -(n_max):(n_max)])

        ### plot bare truncated series
        label=""
        if (koπ_pos==1 && m==0)  label = L"\mathrm{bare \; series} \; n_{max}="*string(n_max) end    
        plot!(plt_vec[m_pos],x_vec_bare,p_x.(x_vec_bare),label=label,color=:grey,linewidth=0.5)
        pade_vec = [[6,6],[5,5]]

        ### x-Padé approximants
        label=""
        for  (pade_pos,pade) in enumerate(pade_vec)
            if koπ_pos==1 && m==1 label="x-Padé "*string(pade) end
            y_vec = get_pade(p_x,pade[1],pade[2]).(x_vec)
            x_vec_show = 1*x_vec
            y_vec_show = 1*y_vec
            for p in 1:length(y_vec)
                if abs(y_vec_show[p])>0.75
                    x_vec_show = 1*x_vec_show[1:p]
                    y_vec_show = 1*y_vec_show[1:p]
                    break
                end
            end   
            plot!(plt_vec[m_pos],x_vec_show,y_vec_show,color=:darkblue,alpha=0.6,linestyle=linestyle_vec[pade_pos],label=label)
        end

        ### u-Padé approximants
        f=0.205
        ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
        p_u = Polynomial(ufromx_mat*coeffs(p_x))
        u_vec = tanh.(f .* x_vec)
        label=""
        for  (pade_pos,pade) in enumerate(pade_vec)
            if koπ_pos==1 && m==2 label="u-Padé "*string(pade)*" (f=$f)" end
            plot!(plt_vec[m_pos],x_vec,get_pade(p_u,pade[1],pade[2]).(u_vec),color=:purple,alpha=0.6,linestyle=linestyle_vec[pade_pos],label=label)
        end
    end

    ### plot QMC data
    for (koπ_pos,koπ) in enumerate(koπ_vec), (m_pos,m) in enumerate(m_vec)
        TGk_m_QMC = [sum([cos(di*koπ*π)*TGiip_m_QMC[x_pos,1+abs(di),1+m] for di in -imax:1:imax]) for x_pos in eachindex(x_vec_QMC)]
        if !(m!=0 && koπ==0) scatter!(plt_vec[m_pos],x_vec_QMC,TGk_m_QMC,color=color_vec[7+koπ_pos],label="",markersize=6,marker=:square) end
        if m==0 annotate!(plt_vec[m_pos],x_vec_QMC[2]+0.1,TGk_m_QMC[2]+0.1-(koπ_pos-1)*0.18,text("k/π="*string(koπ),7,color_vec[7+koπ_pos])) end
    end

    ### put panels together
    xPlots,yPlots=1,3
    plt_final = plot(plt_vec...,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.4*aps_width*yPlots), dpi=600)
    display(plt_final)
    savefig(plt_final,"CaseStudy/Chain_BSb/Chain_S1half_Matsubara_k-space.png")
end