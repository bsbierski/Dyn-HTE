######### Dyn-HTE for chain #########
using Polynomials, HDF5, Measurements
include("plotConventions.jl") 

path_DynHTSE="C:/Users/ruben/Documents/GitHub/Projects/Master/Dyn-HTE/Dyn-HTE/"
include(path_DynHTSE*"Embedding.jl")
include(path_DynHTSE*"LatticeGraphs.jl")
include(path_DynHTSE*"ConvenienceFunctions.jl")

#specify max order
max_order = 12

#LOAD FILES -------------------------
#generate list of graphs
graphs_vec = [load_object(path_DynHTSE*"GraphFiles_chain/graphs_"*string(nn)*".jld2") for nn in 0:max_order];
gG_vec = getGraphsG(graphs_vec);
 
#create vector of all lower order dictionaries
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,max_order+1) 
   
#load dictionaries of all lower orders C_Dict_vec 
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object(path_DynHTSE*"GraphFiles_chain/GraphG_Lists/C_"*string(ord)*".jld2")
end 
#-----------------------------------

### Define Lattice
lattice,LatGraph,center_sites = getLattice_Ball(max_order,"chain");
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.1,fontsize=7,nodeshape=:rect,curves=false))

##### Compute all non-zero G-coefficients c_ii'(m) 
c_iipDyn_mat = get_c_iipDyn_mat_slow(LatGraph,lattice,center_sites,max_order,gG_vec,C_Dict_vec)

cd("/home/bjoern/Dropbox/Desktop/39_KopietzSpinfRG/Num_DynHTE_results/")


###### load QMC data at m=0,1,2 (have x=βJ=[1,2,...,8]; but one can use less) 
x_QMC_vec = collect(1:8)
imax=20
m_vec = [0,1,2]
TGiip_m = zeros((length(x_QMC_vec),imax+1,length(m_vec)))
TGiip_m_err = zeros((length(x_QMC_vec),imax+1,length(m_vec)))

for x_QMC_pos in eachindex(x_QMC_vec) 
    x_QMC=x_QMC_vec[x_QMC_pos]
    fid = h5open("/home/bjoern/Dropbox/Desktop/39_KopietzSpinfRG/Num_DynHTE_results/QMC_Worm/SpinHalfAFMHeisenbergChain/job_BSb_beta"*string(x_QMC)*".out.h5", "r")
    for (m_pos,m) in enumerate(m_vec)
        TGiip_m[x_QMC_pos,:,m_pos] = read(fid["simulation"]["results"]["DensDens_CorrFun_w$m"]["mean"]["value"])[1:imax+1]
        TGiip_m_err[x_QMC_pos,:,m_pos] = read(fid["simulation"]["results"]["DensDens_CorrFun_w$m"]["mean"]["error"])[1:imax+1]
    end
end

###### plot QMC data against Dyn-HTE (r-space)
x_vec = collect(0:0.05:8.1)
x_vec_bare = collect(0:0.05:2)


### PLOT m=0,1,2 Matsubara data against QMC
if false
    di_vec = [0,1,2,4]
    
    ### prepare panels for m=0,1,2
    plt_m0=plot([0,maximum(x_vec)],[0,0],color=:black,lw=0.5,xlims=(0,maximum(x_vec)),ylims=(-0.11,0.26),xformatter=:none,ylabel=L"TG_{ii^\prime}(i\nu_{m=0})", label="",legend=:topright)
    plt_m1=plot([0,maximum(x_vec)],[0,0],color=:black,lw=0.5,xlims=(0,maximum(x_vec)),ylims=(-0.025,0.04),xformatter=:none,ylabel=L"TG_{ii^\prime}({i\nu_{m=1}})", label="",legend=:topleft)
    plt_m2=plot([0,maximum(x_vec)],[0,0],color=:black,lw=0.5,xlims=(0,maximum(x_vec)),ylims=(-0.01,0.017),xlabel=L"x=\beta / J",ylabel=L"TG_{ii^\prime}(i\nu_{m=2})", label="",legend=:topleft)
    plt_vec = [plt_m0,plt_m1,plt_m2]

    ### plot QMC data
    for di in di_vec, (m_pos,m) in enumerate(m_vec)
        scatter!(plt_vec[m_pos],x_QMC_vec,TGiip_m[:,1+di,1+m],color=color_vec[di+1],label="",markersize=5)
        if m==0 annotate!(plt_vec[m_pos],x_QMC_vec[1]+0.5,TGiip_m[3,1+di,1+m],text("i-i'="*string(di),7,color_vec[di+1])) end
    end

    ### plot Dyn-HTE results
    for di in di_vec, (m_pos,m) in enumerate(m_vec)
        
        ### prepare series as polynomial in x
        p = get_TGiip_m_bare(c_iipDyn_mat,m,12)[center_sites[1]+di,1]
        
        ### plot bare truncated series
        label=""
        if (di==0 && m==0)  label = "bare n="*string(max_order) end    
        plot!(plt_vec[m_pos],x_vec_bare,p.(x_vec_bare),label=label,color=:grey)

        ### Padé approximants
        for (pades_pos,pades) in enumerate([[6,6],[5,5]])
            label=""
            if (m_pos==2 && di==0) label="Padé "*string(pades) end
            plot!(plt_vec[m_pos],x_vec,get_pade(p,pades...).(x_vec),linestyle=linestyle_vec[pades_pos],label=label,color=:darkblue,alpha=0.5)
        end
    end

    ### integrated-differential approximants
    for di in di_vec, (m_pos,m) in enumerate(m_vec[1])

        p = get_TGiip_m_bare(c_iipDyn_mat,m,12)[center_sites[1]+di,1]
        
        for (IDA_pos,IDA) in enumerate([[2,4,4]])
            label=""
            if (m_pos==0 && di==0) label="IDA "*string(IDA) end
            plot!(plt_vec[m_pos],x_vec,get_intDiffApprox(p,x_vec,IDA...),linestyle=linestyle_vec[IDA_pos],label=label,color=:magenta,alpha=0.5)
        end
    end



    ### put panels together
    xPlots,yPlots=1,3
    plt_final = plot(plt_vec...,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.42*aps_width*yPlots))
    display(plt_final)
    #savefig(plt_final,"HeisenbergAFMSpinHalfChain_Gr_m_withQMC_n"*string(max_order)*"_tmp.svg")
end




#########  RUBEN STARTED HERE   #############

###### EXTRAPOLATION OF MOMENTS TO LOW TEMPERATURES
x_range = collect(0:0.01:3.5)
k_vec=[0.5*π] #0.8*pi


#plt_m_x = [plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",xlim=(0,3.5),title="r="*string(i-1)) for i=1:6 ]
plt_m_x = [plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",xlim=(0,3.5),title="r="*string(0),legend=:bottomleft),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",xlim=(0,3.5),title="r="*string(1),legend=:bottomleft),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",xlim=(0,3.5),title="r="*string(2),legend=:topright),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",xlim=(0,3.5),title="r="*string(3),legend=:topright),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",xlim=(0,3.5),title="r="*string(4),legend=:topright),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",xlim=(0,3.5),title="r="*string(5),legend=:topright) ]

for (k_pos,k) in enumerate(k_vec)
    c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
    m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)

    c_kDyn_mat1 = get_c_kDyn_mat([(0.8*pi,0)],c_iipDyn_mat,lattice,center_sites)[1]
    m_vec1 = get_moments_from_c_kDyn_mat(c_kDyn_mat1)

    for m_idx=1:length(m_vec)-1
        if m_idx ==1
            plot!(plt_m_x[m_idx],collect(0:0.01:1.9),m_vec[m_idx].(collect(0:0.01:1.9)),label=L"k=0.5 \pi",alpha= 0.5,color = color_vec[1], order=1)
            plot!(plt_m_x[m_idx],collect(0:0.01:1.9),m_vec1[m_idx].(collect(0:0.01:1.9)),label=L"k=0.8 \pi",alpha= 0.5,color = color_vec[2], order=2)
        else 
            plot!(plt_m_x[m_idx],collect(0:0.01:1.9),m_vec[m_idx].(collect(0:0.01:1.9)),label=nothing,alpha= 0.5,color = color_vec[1], order=1)
            plot!(plt_m_x[m_idx],collect(0:0.01:1.9),m_vec1[m_idx].(collect(0:0.01:1.9)),label=nothing,alpha= 0.5,color = color_vec[2], order=2)
        end
        #now pade
        pade_approximant = get_pade(m_vec[m_idx],7-m_idx,7-m_idx)
        plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label="Pade["*string(7-m_idx)*","*string(7-m_idx)*"]",linestyle=linestyle_vec[4],alpha= 0.6,color = color_vec[3])
        #other pades 
        if m_idx<5
            pade_approximant = get_pade(m_vec[m_idx],8-m_idx,6-m_idx)
            plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label="Pade["*string(8-m_idx)*","*string(6-m_idx)*"]",alpha= 1,linestyle=linestyle_vec[2],color = color_vec[4])
            pade_approximant = get_pade(m_vec[m_idx],6-m_idx,8-m_idx)
            plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label="Pade["*string(6-m_idx)*","*string(8-m_idx)*"]",alpha= 1,linestyle=linestyle_vec[3],color = color_vec[5])
        elseif m_idx ==6
            pade_approximant = get_pade(m_vec[m_idx],8-m_idx,6-m_idx)
            plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label="Pade["*string(8-m_idx)*","*string(6-m_idx)*"]",alpha= 1,linestyle=linestyle_vec[2],color = color_vec[4])
            pade_approximant = get_pade(m_vec[m_idx],6-m_idx,8-m_idx)
            plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label="Pade["*string(6-m_idx)*","*string(8-m_idx)*"]",alpha= 1,linestyle=linestyle_vec[3],color = color_vec[5])
        elseif m_idx ==5
            x_range_pade =collect(0:0.01:1.6)
            pade_approximant = get_pade(m_vec[m_idx],8-m_idx,6-m_idx)
            plot!(plt_m_x[m_idx],x_range_pade,pade_approximant.(x_range_pade),label="Pade["*string(8-m_idx)*","*string(6-m_idx)*"]",alpha= 1,linestyle=linestyle_vec[2],color = color_vec[4])
            pade_approximant = get_pade(m_vec[m_idx],6-m_idx,8-m_idx)
            plot!(plt_m_x[m_idx],x_range_pade,pade_approximant.(x_range_pade),label="Pade["*string(6-m_idx)*","*string(8-m_idx)*"]",alpha= 1,linestyle=linestyle_vec[3],color = color_vec[5])
        end
        #now IDA 
        #IDA_approximant = get_intDiffApprox(m_vec[m_idx],x_range,2,4,4)
        #plot!(plt_m_x[m_idx],x_range,IDA_approximant.(x_range),label="IDA[2,4,4]",alpha= 1)
    end

end

k_vec=[0.8*pi]
for (k_pos,k) in enumerate(k_vec)
    c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
    m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)

    for m_idx=1:length(m_vec)-1
        #plot!(plt_m_x[m_idx],collect(0:0.01:1.9),m_vec[m_idx].(collect(0:0.01:1.9)),label=L"k=0.8 \pi",alpha= 0.3)
        #now pade
        pade_approximant = get_pade(m_vec[m_idx],7-m_idx,7-m_idx)
        Plots.plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label=nothing,linestyle=linestyle_vec[4],alpha= 0.6,color = color_vec[3])
        #other pades 
        if m_idx<5
            pade_approximant = get_pade(m_vec[m_idx],8-m_idx,6-m_idx)
            Plots.plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label=nothing,alpha= 1,linestyle=linestyle_vec[2],color = color_vec[4])
            pade_approximant = get_pade(m_vec[m_idx],6-m_idx,8-m_idx)
            Plots.plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label=nothing,alpha= 1,linestyle=linestyle_vec[3],color = color_vec[5])
        elseif m_idx ==6
            pade_approximant = get_pade(m_vec[m_idx],8-m_idx,6-m_idx)
            Plots.plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label=nothing,alpha= 1,linestyle=linestyle_vec[2],color = color_vec[4])
            pade_approximant = get_pade(m_vec[m_idx],6-m_idx,8-m_idx)
            Plots.plot!(plt_m_x[m_idx],x_range,pade_approximant.(x_range),label=nothing,alpha= 1,linestyle=linestyle_vec[3],color = color_vec[5])
        elseif m_idx ==5
            x_range_pade =collect(0:0.01:1.6)
            pade_approximant = get_pade(m_vec[m_idx],8-m_idx,6-m_idx)
            Plots.plot!(plt_m_x[m_idx],x_range_pade,pade_approximant.(x_range_pade),label=nothing,alpha= 1,linestyle=linestyle_vec[2],color = color_vec[4])
            pade_approximant = get_pade(m_vec[m_idx],6-m_idx,8-m_idx)
            Plots.plot!(plt_m_x[m_idx],x_range_pade,pade_approximant.(x_range_pade),label=nothing,alpha= 1,linestyle=linestyle_vec[3],color = color_vec[5])
        end
        #now IDA 
        #IDA_approximant = get_intDiffApprox(m_vec[m_idx],x_range,2,4,4)
        #plot!(plt_m_x[m_idx],x_range,IDA_approximant.(x_range),label="IDA[2,4,4]",alpha= 1)
    end

end

display(plot(plt_m_x...))



###### EXTRAPOLATION OF DELTAS TO INFINITY
x= 0.0
k_vec=[0.5*π]
plt_δ = plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_r",legend=:topleft,title="X="*string(x),xlim=(-0.2,10.2),ylim=(-0.2,8))

ext_vec = [[3,4,true],[3,4,false],[2,4,true],[2,4,false]]

for (k_pos,k) in enumerate(k_vec)
    c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
    m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)[1:6]
    δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec])

    Plots.scatter!(plt_δ,0:(length(δ_vec)-1),δ_vec,color=color_vec[k_pos],label=L"k/\pi="*string(k/π),markersize=7.0,markeralpha=0.6)

    for (ext_pos,ext) in enumerate(ext_vec)
        δ_vec_ext =  extrapolate_δvec(δ_vec,ext[1],ext[2],400,Bool(ext[3]))
        Plots.scatter!(plt_δ,0:length(δ_vec_ext)-1,δ_vec_ext,markersize=3,label="",marker=marker_vec[ext_pos],color=color_vec[k_pos])
    end

end
plt_final = plot(plt_δ)
display(plt_final)

savefig(plt_δ,"HeisenbergAFMSpinHalfChain_moment_extrapolation.svg")


###### S(k,w): few k points, plot δ_vec and show robustness of δ-extrapolation
### pick x=0 for now (later check x>0, use PADE or IDA), get δs and plot them
x=0.0
k_vec=[0.2*π,0.8*π]

w_vec = collect(-3.0:0.01:3.0)
η=0.01
ext_vec = [[2,3,true],[2,3,false],[2,4,true],[2,4,false]]

title=""#"AFM Heisenberg chain, J/T=x=$x"
plt_δ = plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_r",title="X="*string(x),legend=:topleft,xlim=(-0.2,10.2),ylim=(0,8))
plt_S = plot(xlabel=L"w=\omega/J",ylabel=L"JS(k,\omega)",title="X="*string(x))
#plt_m = plot([0],[0],label="",xlabel=L"r",ylabel=L"m_r",xlim=(-0.2,10.2))

for (k_pos,k) in enumerate(k_vec)

    c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
    m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)
    δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec])

    scatter!(plt_δ,0:(length(δ_vec)-1),δ_vec,color=color_vec[k_pos],label=L"k/\pi="*string(k/π),markersize=7.0,markeralpha=0.6)
    
    #scatter!(plt_m,0:(length(δ_vec)-1),[m(x) for m in m_vec],color=color_vec[k_pos],label=L"k/\pi="*string(k/π),markersize=7.0,markeralpha=0.6)

    for (ext_pos,ext) in enumerate(ext_vec)
        δ_vec_ext =  extrapolate_δvec(δ_vec,ext[1],ext[2],2000,Bool(ext[3]))
        scatter!(plt_δ,0:length(δ_vec_ext)-1,δ_vec_ext,markersize=3,label="",marker=marker_vec[ext_pos],color=color_vec[k_pos])
        plot!(plt_S,w_vec,[JS(δ_vec_ext ,x,w,η) for w in w_vec],label="",color=color_vec[k_pos],linestyle=linestyle_vec[ext_pos])
    end
end

xPlots,yPlots=2,1
plt_final = plot(plt_δ,plt_S , layout=(yPlots,xPlots), size=(aps_width*xPlots,0.62*aps_width*yPlots))
display(plt_final)
savefig(plt_final,"HeisenbergAFMSpinHalfChain_Tinf_Skw_various_k.svg")


###### S(k,w) heatmap
using CairoMakie

x = 2.048
k_vec = [(k,0.0) for k in 0.01:0.0039*1.2:(2*π-0.01)]
w_vec = collect(-3.0:0.0314*1.2:3.0)
JSkw_mat = get_JSkw_mat_finitex(x,k_vec,w_vec,0.01,2,3,2000,false,c_iipDyn_mat,lattice,center_sites)


fig = Figure(size=(aps_width*xPlots,0.8*aps_width*yPlots),fontsize=25)
ax=Axis(fig[1,1],limits=(0,2,-3,3),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title="X="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25)
hm=CairoMakie.heatmap!(ax,[k[1]/π for k in k_vec],w_vec,JSkw_mat,colormap=:viridis,colorrange=(0.0,0.25),highclip=:white)
#CairoMakie.Colorbar(fig[:, end+1], hm,size=40)
resize_to_layout!(fig)
display(fig)
#now calculate equal time correlator =1/4 for consitency
equal_time_corr = 0.0
for k_idx in 1:length(k_vec) 
    for w_idx in 1:length(w_vec) 
        equal_time_corr+= JSkw_mat[k_idx,w_idx]
    end 
end
equal_time_corr *= (maximum(w_vec)-minimum(w_vec))/(length(k_vec)*length(w_vec))
println("Equal time correlator: "*string(equal_time_corr))
println("Deviation from 1/4: "*string(round(abs(1/4-equal_time_corr)*4*100,digits=4))*"%")

save("plots/HTSE"*string(x)*"_hq.pdf",fig)


####### S(k,w) heatmap DMRG data
using MAT
beta = 0.128*0#0.128*16
data = matread("SF_beta_"*string(beta)*".mat")

k_DMRG = data["kk"]
w_DMRG = data["om"]
S_DMRG = data["Skw"]

#blur data
using ImageFiltering
sigma = 2.1
kernel = Kernel.gaussian(sigma)
S_DMRG_blurred = imfilter(S_DMRG, kernel)


fig = Figure(size = 2,fontsize=10)
ax=Axis(fig[1,1] ,limits=(0,2,-3,3),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title=L"x=0.0,\:  \sigma=",titlesize=15,xlabelsize=12,ylabelsize=12)

hm=CairoMakie.heatmap!(ax,k_DMRG,w_DMRG,S_DMRG_blurred,colormap=:viridis,colorrange=(0.0,0.25),ylim=(-3,3),highclip=:white)

resize_to_layout!(fig)
display(fig)
#now calculate equal time correlator =1/4 for consitency
equal_time_corr = 0.0
for k_idx in 1:length(k_DMRG) 
    for w_idx in 1:length(w_DMRG) 
        equal_time_corr+= S_DMRG[k_idx,w_idx]
    end 
end
equal_time_corr *= 0.5*(maximum(k_DMRG)-minimum(k_DMRG))*(maximum(w_DMRG)-minimum(w_DMRG))/(length(k_DMRG)*length(w_DMRG))
println("Equal time correlator: "*string(equal_time_corr))
println("Deviation from 1/4: "*string(round(abs(1/4-equal_time_corr)*4*100,digits=4))*"%")

save("plots/DMRG1"*string(beta)*".pdf",fig)


######NOW THE COMPARISON OF Dyn-HTE and DMRG
using MAT
betas = [0.0,1.024,2.048]

title_string =[L"x=0.0,\:  \sigma=0.4%",L"x=1.024,\:  \sigma=0.4%",L"x=2.048,\:  \sigma=0.5%"
,L"x=0.0,\:  \sigma= 1.4%",L"x=1.024,\:  \sigma=0.9%",L"x=2.048,\:  \sigma= 0.7%"] 


fig = CairoMakie.Figure(resolution = (900,500),fontsize=10)

for (i,beta) in enumerate(betas)
    data = matread("SF_beta_"*string(beta)*".mat")

    #DMRG
    k_DMRG = data["kk"]
    w_DMRG = data["om"]
    S_DMRG = data["Skw"]

    #dyn-HTE
    k_vec = [k for k in 0.01:0.0039*1.2:(2*π-0.01)]
    w_vec = collect(-3.0:0.0314*1.2:3.0)

    #blur data
    using ImageFiltering
    sigma = 2.1
    kernel = Kernel.gaussian(sigma)
    S_DMRG_blurred = imfilter(S_DMRG, kernel)




    ax=CairoMakie.Axis(fig[1,i],limits=(0,2,-3,3),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title=title_string[i],titlesize=15,xlabelsize=12,ylabelsize=12)
    #dynHTE spot
    ax_dynHTE=CairoMakie.Axis(fig[2,i],limits=(0,2,-3,3),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title=title_string[i+3],titlesize=15,xlabelsize=12,ylabelsize=12)

    hm=CairoMakie.heatmap!(ax,vec(k_DMRG)[1:2:end],vec(w_DMRG)[1:2:end],S_DMRG_blurred[1:2:end,1:2:end],colormap=:viridis,colorrange=(0.0,0.25),ylim=(-3,3),highclip=:white)
    
    #insert Dyn-HTE ax_dynHTE[i]
    # if i==2
    #     i = 3
    # elseif (i==3)
    #     i=2
    # end
    CairoMakie.heatmap!(ax_dynHTE,vec([k/π for k in k_vec])[1:2:end],vec(w_vec)[1:2:end],S_dyn_HTE[i][1:2:end,1:2:end],colormap=:viridis,colorrange=(0.0,0.25),ylim=(-3,3),highclip=:white)
end



subgrid = GridLayout(fig[1, 4], tellheight = false)
subgrid2 = GridLayout(fig[2, 4], tellheight = false)

Label(subgrid[1, 1], L"JS(\mathbf{k},\omega)",fontsize=14)
Label(subgrid2[1, 1],L"JS(\mathbf{k},\omega)",fontsize=14)

cb1=CairoMakie.Colorbar(subgrid[2, 1],hm,size=20,labelsize = 22)
cb2= CairoMakie.Colorbar(subgrid2[2, 1],hm,size=20,labelsize = 22)



#resize_to_layout!(fig)
#Plots.savefig("testt.pdf")
display(fig)



#now calculate equal time correlator =1/4 for consitency
equal_time_corr = 0.0
for k_idx in 1:length(k_DMRG) 
    for w_idx in 1:length(w_DMRG) 
        equal_time_corr+= S_DMRG[k_idx,w_idx]
    end 
end
equal_time_corr *= 0.5*(maximum(k_DMRG)-minimum(k_DMRG))*(maximum(w_DMRG)-minimum(w_DMRG))/(length(k_DMRG)*length(w_DMRG))
println("Equal time correlator: "*string(equal_time_corr))
println("Deviation from 1/4: "*string(round(abs(1/4-equal_time_corr)*4*100,digits=4))*"%")

save("plots/DMRG1"*string(beta)*".pdf",fig)


############Compare slices of heatmap in DMRG and Dyn-HTE
x=0.128*16
k_vec=[0.2*π,0.8*π]

plt_S = plot(xlabel=L"w=\omega/J",ylabel=L"JS(k,\omega)",title="X="*string(x),titlesize=30)

w_vec = collect(-3.0:0.01:3.0)
η=0.01
ext_vec = [[2,3,false]]

#plt_m = plot([0],[0],label="",xlabel=L"r",ylabel=L"m_r",xlim=(-0.2,10.2))
#load DMRG 
using MAT
data = matread("SF_beta_"*string(x)*".mat")
k_DMRG = data["kk"]
w_DMRG = data["om"]
S_DMRG = data["Skw"]

#blur data
using ImageFiltering
sigma = 2.1
kernel = Kernel.gaussian(sigma)
S_DMRG_blurred = imfilter(S_DMRG, kernel)

for (k_pos,k) in enumerate(k_vec)
    k_slice =  S_DMRG_blurred[round(Int,k/(pi*0.00390625)),:][1:191]
    if k_pos == 1
        Plots.scatter!(plt_S,w_DMRG[1:191],k_slice,color=color_vec[k_pos],label=L"DMRG, \ \  k= 0.2 \pi",alpha=0.1)
    else
        Plots.scatter!(plt_S,w_DMRG[1:191],k_slice,color=color_vec[k_pos],label=L"DMRG,\ \  k= 0.8 \pi",alpha=0.1)
    end
end


for (k_pos,k) in enumerate(k_vec)

    c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
    m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)[1:7]
    m_vec_extrapolated = [get_pade(m_vec[m_idx],7-m_idx,7-m_idx) for m_idx=1:7]
    δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec_extrapolated])

    #scatter!(plt_δ,0:(length(δ_vec)-1),δ_vec,color=color_vec[k_pos],label=L"k/\pi="*string(k/π),markersize=7.0,markeralpha=0.6)
    
    #scatter!(plt_m,0:(length(δ_vec)-1),[m(x) for m in m_vec],color=color_vec[k_pos],label=L"k/\pi="*string(k/π),markersize=7.0,markeralpha=0.6)

    for (ext_pos,ext) in enumerate(ext_vec)
        δ_vec_ext =  extrapolate_δvec(δ_vec,ext[1],ext[2],6000,Bool(ext[3]))
        #scatter!(plt_δ,0:length(δ_vec_ext)-1,δ_vec_ext,markersize=3,label="",marker=marker_vec[ext_pos],color=color_vec[k_pos])
        if ext_pos ==1
            if k_pos ==1
                Plots.plot!(plt_S,w_vec,[JS(δ_vec_ext ,x,w,η) for w in w_vec],color=color_vec[k_pos+2],linestyle=linestyle_vec[ext_pos],alpha=1,label=L"Dyn-HTE, \ \ k= 0.2 \pi")
            else 
                Plots.plot!(plt_S,w_vec,[JS(δ_vec_ext ,x,w,η) for w in w_vec],color=color_vec[k_pos+2],linestyle=linestyle_vec[ext_pos],alpha=1,label=L"Dyn-HTE, \ \  k=0.8 \pi ")
            end
        else 
            Plots.plot!(plt_S,w_vec,[JS(δ_vec_ext ,x,w,η) for w in w_vec],color=color_vec[k_pos+2],linestyle=linestyle_vec[ext_pos],alpha=1,label=nothing)
        end
    end
end



display(plt_S)



#############################################################
#############################################################
#############################################################





###### k-space plots (kind of old)
if false

    function bareSeries_k(L::SimpleGraph,k::Float64,m::Int,max_order::Int)::Polynomial
        """ k in BZ """
        return sum([2*cos(k*r)*bareSeries_r(L,[20,20+r],m,max_order) for r in 0:max_order])-bareSeries_r(L,[20,20],m,max_order)
    end
    
    
    k_vec = [0.0,0.8*π,π]
    ### m=0 plot
    m=0
    plt=plot([0,maximum(x_vec)],[0,0],color=:black,xlims=(0,maximum(x_vec)),ylims=(0.0,0.75),xformatter=:none,ylabel=L"TG_{k}(i\nu_{m=0})", label="",legend=:topleft)
    for k_pos in eachindex(k_vec)
        k = k_vec[k_pos]          
        p = bareSeries_k(L,k,m,max_order)
        pade_vec = [Polynomials.PolyCompat.PadeApproximation.Pade(p,padeNM[1],padeNM[2]) for padeNM in padeNM_vec]
        plot!(plt,x_vec_bare,p.(x_vec_bare),label="bare n="*string(max_order),color=:grey)
        for pade_pos in eachindex(padeNM_vec)
            pade = pade_vec[pade_pos]
            padeNM = padeNM_vec[pade_pos]
            plot!(plt,x_vec,pade.(x_vec),label="",color=:darkblue,linestyle=linestyle_vec[pade_pos])
        end
    end
    for l in eachindex(plt.series_list)
        if l>1+length(padeNM_vec)+1
            plt.series_list[l][:label]=""
        end
    end
    ### QMC data 
    for k_pos in eachindex(k_vec)
        k = k_vec[k_pos] 
        Gk_x_vec = sum([2*cos(k*r)*TGiip_m[:,1+r,1] for r in 0:length(TGiip_m[1,:,1])-1])-TGiip_m[:,1,1]
        scatter!(plt,x_QMC_vec,Gk_x_vec,color=color_vec[end-k_pos],markershape=:rect,label="",markersize=5)
        annotate!(plt,x_QMC_vec[3],Gk_x_vec[2],text("k/π="*string(k/π),7,color_vec[end-k_pos]))
    end
    plt_k_m0 = plt

    ### m=1 plot (avoid k=0 due to total Sz conservation)
    m=1
    plt=plot([0,maximum(x_vec)],[0,0],color=:black,xlims=(0,maximum(x_vec)),ylims=(0.0,0.1),xformatter=:none,ylabel=L"TG_{k}(i\nu_{m=1})", label="",legend=:topleft)
    for k_pos in 2:3
        k = k_vec[k_pos]          
        p = bareSeries_k(L,k,m,max_order)
        pade_vec = [Polynomials.PolyCompat.PadeApproximation.Pade(p,padeNM[1],padeNM[2]) for padeNM in padeNM_vec]
        plot!(plt,x_vec_bare,p.(x_vec_bare),label="",color=:grey)
        for pade_pos in eachindex(padeNM_vec)
            pade = pade_vec[pade_pos]
            padeNM = padeNM_vec[pade_pos]
            plot!(plt,x_vec,pade.(x_vec),label="Padé ["*string(padeNM[1])*","*string(padeNM[2])*"]",color=:darkblue,linestyle=linestyle_vec[pade_pos])
        end

    end
    for l in eachindex(plt.series_list)
        if l>1+length(padeNM_vec)+1
            plt.series_list[l][:label]=""
        end
    end
    ### QMC data 
    for k_pos in 2:3
        k = k_vec[k_pos] 
        Gk_x_vec = sum([2*cos(k*r)*TGiip_m[:,1+r,2] for r in 0:length(TGiip_m[1,:,2])-1])-TGiip_m[:,1,2]
        scatter!(plt,x_QMC_vec,Gk_x_vec,color=color_vec[end-k_pos],markershape=:rect,label="",markersize=5)
    end
    plt_k_m1 = plt

    ### m=2 plot (avoid k=0 due to total Sz conservation)
    m=2
    plt=plot([0,maximum(x_vec)],[0,0],color=:black,xlims=(0,maximum(x_vec)),ylims=(0.0,0.04),xlabel=L"x=\beta / J",ylabel=L"TG_{k}(i\nu_{m=2})", label="",legend=:topleft)
    for k_pos in 2:3
        k = k_vec[k_pos]          
        p = bareSeries_k(L,k,m,max_order)
        pade_vec = [Polynomials.PolyCompat.PadeApproximation.Pade(p,padeNM[1],padeNM[2]) for padeNM in padeNM_vec]
        plot!(plt,x_vec_bare,p.(x_vec_bare),label="",color=:grey)
        for pade_pos in eachindex(padeNM_vec)
            pade = pade_vec[pade_pos]
            padeNM = padeNM_vec[pade_pos]
            plot!(plt,x_vec,pade.(x_vec),label="",color=:darkblue,linestyle=linestyle_vec[pade_pos])
        end
    end
    for l in eachindex(plt.series_list)
        if l>1+length(padeNM_vec)+1
            plt.series_list[l][:label]=""
        end
    end
    ### QMC data 
    for k_pos in 2:3
        k = k_vec[k_pos] 
        Gk_x_vec = sum([2*cos(k*r)*TGiip_m[:,1+r,3] for r in 0:length(TGiip_m[1,:,3])-1])-TGiip_m[:,1,3]
        scatter!(plt,x_QMC_vec,Gk_x_vec,color=color_vec[end-k_pos],markershape=:rect,label="",markersize=5)
    end
    plt_k_m2 = plt    

    ### put panels together
    xPlots,yPlots=1,3
    plt_final = plot(plt_k_m0,plt_k_m1,plt_k_m2,  layout=(yPlots,xPlots), size=(aps_width*xPlots,0.42*aps_width*yPlots))
    display(plt_final)
    savefig(plt_final,"HeisenbergAFMSpinHalfChain_Gk_m_withQMC_n"*string(max_order)*".svg")
end


