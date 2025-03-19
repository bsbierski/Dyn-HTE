######### Dyn-HTE for chain #########
using Polynomials, HDF5, Measurements


path_DynHTSE="C:/Users/ruben/Documents/GitHub/Projects/Master/Dyn-HTE/Dyn-HTE/"
include(path_DynHTSE*"Embedding.jl")
include(path_DynHTSE*"LatticeGraphs.jl")
include(path_DynHTSE*"ConvenienceFunctions.jl")
include(path_DynHTSE*"plotConventions.jl") 

#specify max order
max_order = 11

#LOAD FILES -------------------------
#generate list of graphs
#graphs_vec = [load_object(path_DynHTSE*"GraphFiles_chain/graphs_"*string(nn)*".jld2") for nn in 0:max_order]
gG_vec_unique = getGraphsG(max_order)
 
#create vector of all lower order dictionaries
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,max_order+1) 
   
#load dictionaries of all lower orders C_Dict_vec 
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object(path_DynHTSE*"GraphFiles/GraphG_Lists_S1/C_"*string(ord)*".jld2")
end 
#-----------------------------------

### Define Lattice
LatGraph = complete_graph(4)
reference_site = 1
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.1,fontsize=7,nodeshape=:rect,curves=false))


on_site_correlator = Calculate_Correlator(LatGraph,reference_site,reference_site,max_order,gG_vec,C_Dict_vec)
for ord = 0:max_order
    on_site_correlator[ord+1] *= (-1)^ord
end

m_vec = []
for m_idx = 0:5
    m_summed = Polynomial([0])
    for ord= 2*m_idx:max_order
        if m_idx==0
            m_summed +=  on_site_correlator[ord+1][m_idx+1]*Polynomial(Polynomial(vcat(zeros(ord-2*m_idx),[1])))
        elseif m_idx%2 ==0 
            m_summed -=   on_site_correlator[ord+1][m_idx+1]*Polynomial(Polynomial(vcat(zeros(ord-2*m_idx),[1])))
        else
            m_summed +=  on_site_correlator[ord+1][m_idx+1]*Polynomial(Polynomial(vcat(zeros(ord-2*m_idx),[1])))
        end
    end
    push!(m_vec,m_summed )

end

#moment extrapolation
m_vec_extrapolated_pade = []
for m_idx=1:length(m_vec)
    push!(m_vec_extrapolated_pade, get_pade(m_vec[m_idx],6-m_idx,6-m_idx))
end

plt_m = plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",legend=:topleft,title="moments")
x_vec = 0:0.05:0.5
#plot moments
for m_idx = 1:length(m_vec)
    plot!(plt_m,x_vec,m_vec[m_idx].(x_vec),label="m"*string((m_idx-1)*2))
end
display(plt_m)

x=2.0
δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec_extrapolated_pade])
#plot deltas 
plt_deltas = plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_r",legend=:topleft,title="deltas")
scatter!(plt_deltas,0:length(δ_vec)-1,δ_vec)
display(plt_deltas)

δ_vec_ext = extrapolate_δvec(δ_vec,3,3,1000,true)
w_step_size=0.01
w_vec = collect(-3:w_step_size:3)
JSKw = [JS(δ_vec_ext ,x,w,0.01) for w in w_vec]#[2*pi*(1-exp(-w*x))JS(δ_vec_ext ,x,w,0.01) for w in w_vec]

plot(w_vec,JSKw)

#Calculate moments 






#########  RUBEN STARTED HERE   #############




###### S(k,w): few k points, plot δ_vec and show robustness of δ-extrapolation
### pick x=0 for now (later check x>0, use PADE or IDA), get δs and plot them
x=0.0
k_vec=[0.2*π,π]

w_vec = collect(-3.0:0.01:3.0)
η=0.01
ext_vec = [[5,6,true],[5,6,false],[6,6,true]]

title=""#"AFM Heisenberg chain, J/T=x=$x"
plt_δ = plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_r",title="",legend=:topleft,xlim=(-0.2,10.2),ylim=(0,8))
plt_S = plot(xlabel=L"\omega/J=w",ylabel=L"JS(k,\omega)",title="")
#plt_m = plot([0],[0],label="",xlabel=L"r",ylabel=L"m_r",xlim=(-0.2,10.2))

for (k_pos,k) in enumerate(k_vec)

    c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
    m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)
    δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec])

    Plots.scatter!(plt_δ,0:(length(δ_vec)-1),δ_vec,color=color_vec[k_pos],label=L"k/\pi="*string(k/π),markersize=7.0,markeralpha=0.6)
    
    #scatter!(plt_m,0:(length(δ_vec)-1),[m(x) for m in m_vec],color=color_vec[k_pos],label=L"k/\pi="*string(k/π),markersize=7.0,markeralpha=0.6)

    for (ext_pos,ext) in enumerate(ext_vec)
        δ_vec_ext =  extrapolate_δvec(δ_vec,ext[1],ext[2],2000,Bool(ext[3]))
        Plots.scatter!(plt_δ,7:length(δ_vec_ext)-1,δ_vec_ext[8:end],markersize=3,label="",marker=marker_vec[ext_pos],color=color_vec[k_pos])
        Plots.plot!(plt_S,w_vec,[JS(δ_vec_ext ,x,w,η) for w in w_vec],label="",color=color_vec[k_pos],linestyle=linestyle_vec[ext_pos])
        Plots.plot!(plt_S,w_vec,[JS(δ_vec ,x,w,η) for w in w_vec],label="",color=color_vec[k_pos],linestyle=linestyle_vec[1],alpha=0.13)

    end
end

addABC(plt_δ,"a)")
addABC(plt_S,"b)")

xPlots,yPlots=1,2
plt_final = plot(plt_δ,plt_S , layout=(yPlots,xPlots), size=(aps_width*xPlots,0.32*aps_width*yPlots))
display(plt_final)
savefig(plt_final,"chain_Tinf_Skw_various.pdf")






#CHECK WHICH f VALUE IS BEST 
x=2.0
f_values = 0.01:0.01:1

deviation_average = zeros(length(f_values))
for x =1:0.2:4
for k_spacing in [22.2]#[21.9,21.91,21.92,21.93,21.94,21.95,21.96,21.97,21.98,21.99,22.0,22.01,22.02,22.03,22.04,22.05]
k_vec=  0.015:0.0039*k_spacing:(2*π-0.01)

r_max = 3

deviation_arr = zeros(length(f_values))

for (f_idx,f) in enumerate(f_values)
    println(f_idx)
    summ = 0.0
    
    for (k_pos,k) in enumerate(k_vec)
        c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
        m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)[1:6]
        m_vec_times_x =[m_vec[i]*Polynomial([0,1]) for i=1:length(m_vec)]
        m_vec_extrapolated_pade1 = []
        m_vec_extrapolated_pade2 = []


        for m_idx =1:r_max+1
            p_u = Polynomial(subs_matrix_list[f_idx][m_idx]*coeffs(m_vec_times_x[m_idx]))
            push!(m_vec_extrapolated_pade1, get_pade(p_u,8-m_idx,7-m_idx))
            push!(m_vec_extrapolated_pade2, get_pade(p_u,7-m_idx,8-m_idx))

            #sum+= abs(get_pade(p_u,8-m_idx,7-m_idx)(tanh(f*x))/x-get_pade(p_u,7-m_idx,8-m_idx)(tanh(f*x))/x)
        end
        δ_vec1,r_vec = fromMomentsToδ([m(tanh(f*x))/x for m in m_vec_extrapolated_pade1])
        δ_vec2,r_vec = fromMomentsToδ([m(tanh(f*x))/x for m in m_vec_extrapolated_pade2])
        # ###Now extrapolte deltas
        # δ_vec_ext1 = extrapolate_δvec(δ_vec1,3,3,1000,false)
        # δ_vec_ext2 = extrapolate_δvec(δ_vec2,3,3,1000,false)

        # if δ_vec_ext1[end] <0 || δ_vec_ext2[end] <0
        #     summ=-(r_max*length(k_vec))
        #     break
        # else
        #     summ += sum( [ abs( JS(δ_vec_ext1 ,x,w,0.01)-JS(δ_vec_ext2 ,x,w,0.01) ) for w in collect(-3.0:0.0314*2.4:3.0) ] )
            
        # end

        for i=1:r_max+1
            summ += abs(δ_vec1[i]-δ_vec2[i])^2
        end


    end
    deviation_arr[f_idx] = summ/(r_max*length(k_vec))

end

deviation_average+= 1/20*deviation_arr
end
end

Plots.plot(f_values,deviation_average,ylims=(0.0001,100000),yscale=:log10)
display(Plots.scatter!(f_values,deviation_average,ylims=(0.0001,100000),yscale=:log10))
println("f min is :"*string(0.01+(argmin(deviation_average)-1)*0.01))





#PLOT DELTA PARAMETERS FOR DIFFERENT TEMPERATURES
###### delta(beta)

f= 0.31

plot_lst = []
for (k_idx,k) in enumerate([ 0.2*pi,pi])

betas = 0:0.01:4.5
#k= 2.350798079#0.699*pi



#calculate delta for respective x
c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)

# #BASIC PADE 
# m_vec_extrapolated_pade_basic = []
# for m_idx=1:length(m_vec)-2
#     push!(m_vec_extrapolated_pade_basic, get_pade(m_vec[m_idx],7-m_idx,7-m_idx))
# end

#SUBSTITUTION 
     #0.76
pade_orders = [(8,7),(7,8),(7,6)]
m_vec_times_x = [m_vec[i]*Polynomial([0,1]) for i=1:length(m_vec)]
m_vec_times_x_normalized = [m_vec[i]*Polynomial([0,1])/m_vec[i](0) for i=1:length(m_vec)]
m_vec_extrapolated_pade = [[] for i=1:length(pade_orders)]
for (idx,pade_order) in enumerate(pade_orders)
    for m_idx=1:4
        substitution_matrix = subs_matrix_list[round(Int,(f-0.25)/0.01+1)][m_idx]#get_LinearTrafoToCoeffs_u(15-2*m_idx,f)
        p_u = Polynomial(substitution_matrix*coeffs(m_vec_times_x[m_idx]))
        push!(m_vec_extrapolated_pade[idx], get_pade(p_u,pade_order[1]-m_idx,pade_order[2]-m_idx))
    end
end

delta_T = [[[],[],[],[],[],[]]  for i=1:length(pade_orders)]

for (idx,pade_order) in enumerate(pade_orders)
    for x in [0.5,1.0,2.0,4.0]
        #PADE
        #δ_vec_basic,r_vec_basic = fromMomentsToδ([m(x) for m in m_vec_extrapolated_pade_basic])
        #SUBS
        δ_vec,r_vec = fromMomentsToδ([m(tanh(f*x))/x for m in m_vec_extrapolated_pade[idx]])

        #deltas 
        push!(delta_T[idx][1],δ_vec[1])
        push!(delta_T[idx][2],δ_vec[2])
        push!(delta_T[idx][3],δ_vec[3])
        push!(delta_T[idx][4],δ_vec[4])
        #push!(delta_T[idx][5],δ_vec[5])

    end
end


# #BASIC PADE
# delta_T_basic = [[],[],[],[],[],[]]
# for x in betas
#     #PADE
#     δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec_extrapolated_pade_basic])

#     #deltas 
#     push!(delta_T_basic[1],δ_vec[1])
#     push!(delta_T_basic[2],δ_vec[2])
#     push!(delta_T_basic[3],δ_vec[3])
#     push!(delta_T_basic[4],δ_vec[4])
#     #push!(delta_T_basic[5],δ_vec[5])

# end

if k_idx ==1
    plt_δ = plot([0],[0],label="",xlabel=L"r",ylabel=[L"\delta_r",""][k_idx],legend=:topleft,xlim=(-0.2,3.1),ylim=[(-0.1,6.5),(-0.1,4.1),(-0.1,3.7)][k_idx])
    plt_m = plot([0],[0],label="",xlabel=L"x",ylabel=[L"x \cdot m_{\mathbf{k},r}(x)/m_{\mathbf{k},r}(0)",""][k_idx],title=[L"k=0.2 \pi",L"k=\pi"][k_idx],legend=:topleft,xlim=(-0.2,4.5),ylim=[(-0.05,5.1),(-0.05,5.1),(-0.05,2.5)][k_idx])
else 
    plt_δ = plot([0],[0],label="",xlabel=L"r",ylabel=[L"\delta_r",""][k_idx],legend=:topleft,xlim=(-0.2,3.1),ylim=[(-0.1,6.5),(-0.1,6.5),(-0.1,3.7)][k_idx] ,yticks=([ 0,1,2,3,4,5,6], ["", "", "","", "", "",""]))
    plt_m = plot([0],[0],label="",xlabel=L"x",ylabel=[L"x \cdot m_{\mathbf{k},r}(x)/m_{\mathbf{k},r}(0)",""][k_idx],title=[L"k=0.2 \pi",L"k=\pi"][k_idx],legend=:topleft,xlim=(-0.2,4.5),ylim=[(-0.05,5.1),(-0.05,5.1),(-0.05,2.5)][k_idx] ,yticks=([ 0, 1,2,3,4,5], ["", "", "", "", "", ""]))
end

if k_idx ==1 
    Plots.scatter!(plt_δ,[0],[0],label="x=0.5",color = thermalCol4_vec[1])
    Plots.scatter!(plt_δ,[0],[0],label="x=1.0",color = thermalCol4_vec[2])
    Plots.scatter!(plt_δ,[0],[0],label="x=2.0",color = thermalCol4_vec[3])
    Plots.scatter!(plt_δ,[0],[0],label="x=4.0",color = thermalCol4_vec[4])

    Plots.plot!(plt_m,[0],[0],label="x bare",color = "grey",linestyle = linestyle_vec[1],linewidth=0.4)
    Plots.plot!(plt_m,[0],[0],label="u Padé [7-r,6-r]",color = "grey",linestyle = linestyle_vec[2],alpha =0.5)
    Plots.plot!(plt_m,[0],[0],label="u Padé [6-r,7-r]",color = "grey",linestyle = linestyle_vec[3])
    Plots.plot!(plt_m,[0],[0],label="u Padé [6-r,5-r]",color = "grey",linestyle = linestyle_vec[4])
    #Plots.plot!(plt_m,[0],[0],label="x Padé [6-r,6-r]",color = "grey",linestyle = linestyle_vec[4])


    Plots.plot!(plt_m,[0],[0],label="r=0",color = color_vec[1],linestyle = linestyle_vec[1])
    Plots.plot!(plt_m,[0],[0],label="r=1",color = color_vec[2],linestyle = linestyle_vec[1])
    Plots.plot!(plt_m,[0],[0],label="r=2",color = color_vec[3],linestyle = linestyle_vec[1])
    Plots.plot!(plt_m,[0],[0],label="r=3",color = color_vec[4],linestyle = linestyle_vec[1])
end

#BARE SERIES
for i=1:4
    Plots.plot!(plt_m,betas[1:180],m_vec_times_x_normalized[i].(betas[1:180]),label = nothing,alpha= 0.7,color = color_vec[i],linestyle = linestyle_vec[1],linewidth=0.5)
end

#U PADE
for (idx,pade_order) in enumerate(pade_orders)
    for i=1:4
        #Plots.plot!(plt_δ,betas,delta_T[idx][i],label = nothing,alpha= [0.5,1][idx],color = color_vec[i],linestyle = linestyle_vec[idx+1])
        Plots.plot!(plt_m,betas,m_vec_extrapolated_pade[idx][i].(tanh.(f.*betas))/(m_vec_extrapolated_pade[idx][i].(tanh.(f.*0.001))/0.001),label = nothing,alpha= [0.5,1,1][idx],color = color_vec[i],linestyle = linestyle_vec[idx+1])
    end
end
# #BASIC PADE
# for i=1:4
#     #Plots.plot!(plt_δ,betas,delta_T_basic[i],label=nothing,color = color_vec[i],linestyle = linestyle_vec[4])
#     Plots.plot!(plt_m,betas,m_vec_extrapolated_pade_basic[i].(betas).*betas/(m_vec_extrapolated_pade_basic[i].(0.0)),color = color_vec[i],linestyle = linestyle_vec[4],label = nothing)
# end

#DELTAS 
i=0
for beta_idx in 1:4
    i+=1
    Plots.scatter!(plt_δ,0:3,[delta_T[1][1][beta_idx],delta_T[1][2][beta_idx],delta_T[1][3][beta_idx],delta_T[1][4][beta_idx]],label = nothing,color = thermalCol4_vec[i],markersize=7)
end

push!(plot_lst, plt_δ)
push!(plot_lst, plt_m)

end
xPlots,yPlots=2,2
# display(Plots.plot(plot_lst[2],plot_lst[1],  layout=(yPlots,xPlots), size=(aps_width*2,0.60*aps_width*yPlots)))
display(Plots.plot(plot_lst[2],plot_lst[4],plot_lst[1],plot_lst[3],  layout=(yPlots,xPlots), size=(aps_width*2,0.60*aps_width*yPlots)))

savefig("chain_moments_and_deltas.pdf")






#NOW SPIN STRUCTURE FACTOR HEATMAPS
###### S(k,w) heatmap (Dyn-HTE)
using CairoMakie

f=0.48#0.31

x = 4.0#0.128*32 
k_step_size = 0.0039*22.2
w_step_size = 0.05
k_vec = [(k,0.0) for k in 0.015:k_step_size:(2*π-0.01)]#[(k,0.0) for k in 0.01:0.0039*2.4:(2*π-0.01)]
w_vec = collect(-3:w_step_size:3)# for reliable sigma: collect(-5:0.05*0.1:5)
JSkw_mat = get_JSkw_mat_finitex(f,"pade_subs",x,k_vec,w_vec,0.01,3,3,1000,false,c_iipDyn_mat,lattice,center_sites)


fig = Figure(fontsize=10)
ax=Axis(fig[1,1],xlabel=L"k/ \pi",ylabel=L"\omega/J=w",limits=(0,2,-3,3),title="Dyn-HTE x="*string(x)*" f="*string(f),titlesize=15,xlabelsize=12,ylabelsize=12)
hm=CairoMakie.heatmap!(ax,[k[1]/π for k in k_vec],w_vec,JSkw_mat,colormap=:viridis,colorrange=(0.0,0.45),highclip=:white)
# CairoMakie.Colorbar(fig[:, end+1], hm,size=10)
grid = fig[1,2]= GridLayout()
subgrid = GridLayout(grid[1, 1], tellheight = false)
Label(subgrid[1, 1], L"JS(k,\omega)",fontsize=15)
cb1=CairoMakie.Colorbar(subgrid[2, 1],hm,size=11,labelsize = 10) 
display(fig)
# save("chain_x_4.pdf",fig)




#now calculate equal time correlator =1/4 for consitency
equal_time_corr = 0.0
for k_idx in 1:length(k_vec) 
    for w_idx in 1:length(w_vec) 
        equal_time_corr+= JSkw_mat[k_idx,w_idx]
    end 
end
equal_time_corr *= w_step_size*k_step_size/(2*pi)
println("Equal time correlator: "*string(equal_time_corr))
println("Deviation from 1/4: "*string(round(abs(1/4-equal_time_corr)*4*100,digits=4))*"%")

#save("plots/HTSE"*string(x)*"_hq.pdf",fig)


####### S(k,w) heatmap DMRG data
using MAT
beta = 16#0.128*16
data = matread("SF_beta_"*string(beta)*".mat")

k_DMRG = data["kk"]
w_DMRG = data["om"]
S_DMRG = data["Skw"]

#blur data
using ImageFiltering
sigma = 2.1
kernel = Kernel.gaussian(sigma)
S_DMRG_blurred = imfilter(S_DMRG, kernel)


fig = Figure(fontsize=10)
ax=Axis(fig[1,1] ,limits=(0,2,-3,3),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title="DMRG x"*string(x),titlesize=15,xlabelsize=12,ylabelsize=12)

hm=CairoMakie.heatmap!(ax,[k_DMRG[i] for i=1:length(k_DMRG)],[w_DMRG[i] for i=1:length(w_DMRG)],S_DMRG_blurred,colormap=:viridis,colorrange=(0.0,0.4),highclip=:white)

resize_to_layout!(fig)
display(fig)
#now calculate equal time correlator =1/4 for consitency
equal_time_corr = 0.0
for k_idx in 1:length(k_DMRG) 
    for w_idx in 1:length(w_DMRG) 
        equal_time_corr+= S_DMRG[k_idx,w_idx]
    end 
end
equal_time_corr *= 0.5*0.00390625*0.03138610973165286
println("Equal time correlator: "*string(equal_time_corr))
println("Deviation from 1/4: "*string(round(abs(1/4-equal_time_corr)*4*100,digits=4))*"%")

save("plots/DMRG1"*string(beta)*".pdf",fig)




######NOW THE COMPARISON OF Dyn-HTE and DMRG
S_dyn_HTE= load_object("DSF_with_DynHTE_differnt_temperatures.jdl2")
using MAT
betas = [0,2,4]

title_string =["","","","","","","","",""] 

xPlots,yPlots=3,3

fig = CairoMakie.Figure(layout=(yPlots,xPlots), size=(aps_width*2,0.60*aps_width*2),fontsize=7)

grid = fig[1,1]= GridLayout()

plotsDMRG = []
plotsDynHTE = []
plotsKslices = []
hm= 0
for (i,beta) in enumerate(betas)
    data = matread("SF_beta_"*string(beta)*".mat")

    #DMRG
    k_DMRG = data["kk"]
    w_DMRG = data["om"]
    S_DMRG = data["Skw"]

    #dyn-HTE
    k_vec = [k for k in 0.015:0.0039*22.2:(2*π-0.01)]
    w_step_size = 0.05
    w_vec = collect(-3:w_step_size:3)
    

    #blur data
    using ImageFiltering
    sigma = 2.1
    kernel = Kernel.gaussian(sigma)
    S_DMRG_blurred = imfilter(S_DMRG, kernel)




    ax=CairoMakie.Axis(grid[2,i],limits=(0,2,-2.5,2.5),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title=title_string[i],titlesize=10,xlabelsize=8,ylabelsize=8)
    
    

    #DMRG
    ax_dynHTE=CairoMakie.Axis(grid[1,i],limits=(0,2,-2.5,2.5),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title=title_string[i+3],titlesize=10,xlabelsize=8,ylabelsize=8)

    hm=CairoMakie.heatmap!(ax,vec(k_DMRG)[1:3:end],vec(w_DMRG)[1:2:end],S_DMRG_blurred[1:3:end,1:2:end],colormap=:viridis,colorrange=(0.0,0.45),highclip=:white)
    
    push!(plotsDMRG, ax)
    #dynHTE spot
    CairoMakie.heatmap!(ax_dynHTE,vec([k/π for k in k_vec])[1:end],vec(w_vec)[1:end],S_dyn_HTE[i][1:end,1:end],colormap=:viridis,colorrange=(0.0,0.45),highclip=:white)
    push!(plotsDynHTE, ax_dynHTE)

    ###annotate
    CairoMakie.text!(ax_dynHTE, 0.03, 2.42, text=[L"x=0",L"x=2",L"x=4"][i], fontsize=13, color=:white, align=(:left, :top) )
    CairoMakie.text!(ax, 0.03, 2.42, text=[L"x=0",L"x=2",L"x=4"][i], fontsize=13, color=:white, align=(:left, :top) )
    CairoMakie.text!(ax_dynHTE, 0.03, -1.8, text="Dyn-HTE", fontsize=12, color=:white, align=(:left, :top) )
    CairoMakie.text!(ax, 0.03, -1.8, text="DMRG", fontsize=12, color=:white, align=(:left, :top) )
    CairoMakie.text!(ax_dynHTE, 1.29, -1.81, text=[L"\Sigma=0.251",L"\Sigma=0.251",L"\Sigma=0.248"][i], fontsize=13, color=:white, align=(:left, :top) )
    CairoMakie.text!(ax, 1.29, -1.81, text=L"\Sigma=0.250", fontsize=13, color=:white, align=(:left, :top) )

    #k slices 
    k_slice_values=[(0.2*π,0.0),(1*π,0.0)]
    

    if i==1
        ax_k_sliced =CairoMakie.Axis(grid[3,i],xlabel=L"\omega/J=w",ylabel=L"JS(k,\omega)",title=title_string[i+6],titlesize=10,xlabelsize=8,ylabelsize=8 ,limits=(-3,3,-0.01,0.27))

        for (k_pos,k) in enumerate([k_slice_values[i][1] for i=1:length(k_slice_values)])
            CairoMakie.scatter!(ax_k_sliced,w_DMRG[1:2:191],S_DMRG_blurred[round(Integer,k/(pi*0.00390625)),:][1:2:191],color = color_vec[k_pos],markersize =5, alpha =0.6,label=["DMRG \n (k=0.2π)","DMRG \n (k=π)"][k_pos])
            CairoMakie.lines!(ax_k_sliced,w_vec,S_dyn_HTE[i][round(Int,(k-0.015)/(0.0039*22.2))+1,:],color = color_vec[k_pos],label=["Dyn-HTE \n (k=0.2π)","Dyn-HTE \n (k=π)"][k_pos] )
        end
        subgrid = GridLayout(grid[3, 4], tellheight = false)
        legend =Legend(subgrid[1,2],ax_k_sliced, position  = (0, 1),   tellwidth = false, backgroundcolor = :white, framevisible = false, labelsize = 7  ,anchor = :TopRight          )

        CairoMakie.text!(ax_k_sliced,-2.84, 0.237, text=L"x=0", fontsize=10, color=:black)
    else
        ax_k_sliced =CairoMakie.Axis(grid[3,i],xlabel=L"\omega/J=w",title=title_string[i+6],titlesize=10,xlabelsize=8,ylabelsize=8 ,limits=[(-3,3,-0.012,0.31),(-3,3,-0.03,0.64)][i-1])

        for (k_pos,k) in enumerate([k_slice_values[i][1] for i=1:length(k_slice_values)])
            CairoMakie.scatter!(ax_k_sliced,w_DMRG[1:2:191],S_DMRG_blurred[round(Integer,k/(pi*0.00390625)),:][1:2:191],color = color_vec[k_pos],markersize =5, alpha =0.6)
            CairoMakie.lines!(ax_k_sliced,w_vec,S_dyn_HTE[i][round(Int,(k-0.015)/(0.0039*22.2))+1,:],color = color_vec[k_pos])
        end

        CairoMakie.text!(ax_k_sliced, -2.84, [0.269,0.563][i-1], text=[L"x=2",L"x=4"][i-1], fontsize=10, color=:black )
    end

    ax_k_sliced.xgridvisible =false
    ax_k_sliced.ygridvisible =false
    push!(plotsKslices, ax_k_sliced)
end 




subgrid = GridLayout(grid[1, 4], tellheight = false)
subgrid2 = GridLayout(grid[2, 4], tellheight = false)

Label(subgrid[1, 1], L"JS(k,\omega)",fontsize=10)
Label(subgrid2[1, 1],L"JS(k,\omega)",fontsize=10)


cb1=CairoMakie.Colorbar(subgrid[2, 1],hm,size=11,labelsize = 10) 
cb2= CairoMakie.Colorbar(subgrid2[2, 1],hm,size=11,labelsize = 10)

hideydecorations!(plotsDMRG[2])
hideydecorations!(plotsDMRG[3])
hidedecorations!(plotsDynHTE[2], grid = false)
hidedecorations!(plotsDynHTE[3], grid = false)
hidexdecorations!(plotsDynHTE[1])
hidedecorations!(plotsDynHTE[3], grid = false)
# hideydecorations!(plotsKslices[2])
# hideydecorations!(plotsKslices[3])

colgap!(subgrid, 0)
rowgap!(subgrid, 0)

colgap!(subgrid2, 0)
rowgap!(subgrid2, 0)

colgap!(grid, 7)
rowgap!(grid, 8)


resize_to_layout!(fig)
# Plots.savefig("dynHTE_dmrg_test.pdf")
display(fig)
save("chain_dynHTE_DMRG_comparison.pdf",fig)


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
betas = [2.0]
moments = [(5,6),(3,4),(2,3)]
plot_lst = []
for (x_index,x) in enumerate(betas)
    k_slice_values=[(0.1*π,0.0),(1*π,0.0),(0.55*π,0.0)]
    k_vec = [(k,0.0) for k in 0.01:0.0039*2.2:(2*π-0.01)]
    if x_index !=3
        plt_S = plot(xlabel=L"w=\omega/J",ylabel=L"JS(k,\omega)",title=L"x="*string(x),titlefont=13,size=(400,300),legend=:topleft,fontsize=14,ylims=(-0.05,0.6))
    else
        plt_S = plot(xlabel=L"w=\omega/J",ylabel=L"JS(k,\omega)",title=L"x="*string(x),titlefont=13,size=(400,600),legend=:topleft,fontsize=14)
    end
    w_vec = collect(-3.0:0.0314*2.4:3.0)
    η=0.01
    ext_vec = [[2,3,false]]

    #plt_m = plot([0],[0],label="",xlabel=L"r",ylabel=L"m_r",xlim=(-0.2,10.2))
    #load DMRG 
    using MAT
    data = matread("SF_beta_"*string(2)*".mat")
    k_DMRG = data["kk"]
    w_DMRG = data["om"]
    S_DMRG = data["Skw"]

    #blur data
    using ImageFiltering
    sigma = 2.1
    kernel = Kernel.gaussian(sigma)
    S_DMRG_blurred = imfilter(S_DMRG, kernel)
    
    for (k_pos,k) in enumerate([k_slice_values[i][1] for i=1:length(k_slice_values)])
        k_slice =  S_DMRG_blurred[round(Integer,k/(pi*0.00390625)),:][1:191]
        if k_pos == 1
            Plots.scatter!(plt_S,w_DMRG[1:2:191],k_slice[1:2:end],color=color_vec[k_pos],label=[L"DMRG, \ \  k= 0.1 \pi",nothing,nothing][x_index],alpha=0.1)
        else
            Plots.scatter!(plt_S,w_DMRG[1:2:191],k_slice[1:2:end],color=color_vec[k_pos],label=[L"DMRG,\ \  k=  \pi",nothing,nothing][x_index],alpha=0.1)
        end
    end


    JSkw_mat_pade =  JSkw_mat#get_JSkw_mat_finitex("pade",x,k_slice_values,w_vec,0.01,moments[x_index][1],moments[x_index][2],3000,false,c_iipDyn_mat,lattice,center_sites)

    #JSkw_mat_ida = get_JSkw_mat_finitex("ida",x,k_slice_values,w_vec,0.01,2,3,2000,false,c_iipDyn_mat,lattice,center_sites)

    #JSkw_mat_directdelta = get_JSkw_mat_finitex("directdelta",x,k_slice_values,w_vec,0.01,2,3,2000,false,c_iipDyn_mat,lattice,center_sites)


    count=1
    for (k_pos,k) in k_slice_values

        if count ==1
            Plots.plot!(plt_S,w_vec,JSkw_mat_pade[round(Int64,(k_pos-0.015)/(0.0039*31.4)+1),:],color=color_vec[3],linestyle=linestyle_vec[1],alpha=1,label=[L"Dyn-HTE, \ \ k= 0.1 \pi",nothing,nothing][x_index],fontsize=15)
            #Plots.plot!(plt_S,w_vec,JSkw_mat_ida[k_pos,:],color=color_vec[4],linestyle=linestyle_vec[2],alpha=1,label=L"Dyn-HTE (IDA), \ \ k= 0.2 \pi")
            #Plots.plot!(plt_S,w_vec,JSkw_mat_directdelta[k_pos,:],color=color_vec[5],linestyle=linestyle_vec[3],alpha=1,label=L"Dyn-HTE (directly \delta), \ \ k= 0.2 \pi")
        elseif count ==2 
            Plots.plot!(plt_S,w_vec,JSkw_mat_pade[round(Int64,(k_pos-0.015)/(0.0039*31.4)+1),:],color=color_vec[4],linestyle=linestyle_vec[1],alpha=1,label=[L"Dyn-HTE, \ \ k=  \pi",nothing,nothing][x_index],fontsize=15)
            #Plots.plot!(plt_S,w_vec,JSkw_mat_ida[k_pos,:],color=color_vec[4],linestyle=linestyle_vec[2],alpha=1,label=nothing)
            #Plots.plot!(plt_S,w_vec,JSkw_mat_directdelta[k_pos,:],color=color_vec[5],linestyle=linestyle_vec[3],alpha=1,label=nothing)
        else 
            Plots.plot!(plt_S,w_vec,JSkw_mat_pade[round(Int64,(k_pos-0.015)/(0.0039*31.4)+1),:],color=color_vec[5],linestyle=linestyle_vec[1],alpha=1,label=[L"Dyn-HTE, \ \ k= 0.5 \pi",nothing,nothing][x_index],fontsize=15)
        end
        count+=1

    end

    push!(plot_lst,plt_S)
end

plot_total = plot(plot_lst...,layout = (3,1))

display(plot_total) 
#savefig(plot_total, "plots/k_slice_comparison.pdf")




#######Compare delta extraploation Schemes
x_range_raw = collect(0:0.01:1.5)
x_range_moment_pade =collect(0:0.01:2)
k_vec=[0.1*pi]#[0.5*π,0.8*pi] #0.8*pi

bare_range_vec=[1.7,1.5,1.5,1.4,1.2,0.7]
pade_on_moments_range_vec=[3.5,3.5,3.5,3,2,1.1]
ida_on_moments_range_vec=[3.5,3.5,3.5,3.5,3,2.1]
ida_on_deltas_range_vec=[3.5,3.5,3.5,3.5,3,2.1]

#import DMRG data and calculate deltas 
using MAT
betas = [0.0,0.128*8,0.128*16]
DMRG_deltas = [[zeros(3) for i=1:6] for j=1:length(k_vec) ]
for (beta_idx,beta) in enumerate(betas)
    data = matread("SF_beta_"*string(beta)*".mat")
    k_DMRG = data["kk"]
    w_DMRG = data["om"]
    S_DMRG = data["Skw"]
    #blur data
    using ImageFiltering
    sigma = 2#2.1
    kernel = Kernel.gaussian(sigma)
    S_DMRG_blurred = imfilter(S_DMRG, kernel)

    #calculate moments for specific k
    moments = zeros(7)

    for (k_pos,k) in enumerate(k_vec)
        for m_idx=1:6
            res = 0.0
            for (w_idx,w) in enumerate(w_DMRG)
                if beta*w!=0
                    res += (1-exp(-beta*w))/(beta*w)*(S_DMRG_blurred[round(Integer,k/(pi*0.00390625)),w_idx])*w^((m_idx-1)*2)
                else
                    res += (S_DMRG_blurred[round(Integer,k/(pi*0.00390625)),w_idx])*w^((m_idx-1)*2)
                end
            end
            moments[m_idx] = res*9/length(w_DMRG)
        end

        δ_vec,r_vec = fromMomentsToδ(moments)

        for m_idx=1:6
            DMRG_deltas[k_pos][m_idx][beta_idx] = δ_vec[m_idx]
        end

    end


end




#plt_m_x = [plot([0],[0],label="",xlabel=L"x",ylabel=L"m_r",xlim=(0,3.5),title="r="*string(i-1)) for i=1:6 ]
plt_m_x = [plot([0],[0],label="",xlabel=L"x",ylabel=L"\delta_r",xlim=(-0.03,3.5),title="r="*string(0),legend=:bottomleft,size=(600,600)),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"\delta_r",xlim=(-0.03,3.5),title="r="*string(1),legend=:bottomleft),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"\delta_r",xlim=(-0.03,3.5),title="r="*string(2),legend=:bottomleft),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"\delta_r",xlim=(-0.03,3.5),title="r="*string(3),legend=:bottomleft),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"\delta_r",xlim=(-0.03,3.5),title="r="*string(4),legend=:bottomleft),
    plot([0],[0],label="",xlabel=L"x",ylabel=L"\delta_r",xlim=(-0.03,3.5),title="r="*string(5),legend=:topright) ]

for (k_pos,k) in enumerate(k_vec)
    c_kDyn_mat = get_c_kDyn_mat([(k,0)],c_iipDyn_mat,lattice,center_sites)[1]
    m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)

    δ_vec_raw = fromMomentsToδ(m_vec)

    #extrapolate moments with pade 
    m_vec_extrapolated_pade = []
    for m_idx=1:length(m_vec)
        push!(m_vec_extrapolated_pade, get_pade(m_vec[m_idx],7-m_idx,7-m_idx))
    end

    #extrapolate moments with IDA 
    m_vec_extrapolated_ida = []
    for m_idx=1:length(m_vec)
        for idx=0:length(m_vec[m_idx])-1
            if abs(m_vec[m_idx][idx])< 0.00000000000000001
                m_vec[m_idx][idx] =0
            end
        end
        if m_idx <5
            push!(m_vec_extrapolated_ida, get_intDiffApprox(m_vec[m_idx],collect(0:0.01:ida_on_moments_range_vec[m_idx]),2,5-m_idx,5-m_idx))
        end
    end



    #now iterate over deltas
    for m_idx=1:length(m_vec)-1

        #plot the raw delta
        x_range =collect(0:0.01:bare_range_vec[m_idx])
        plot!(plt_m_x[m_idx],x_range,δ_vec_raw[m_idx].(x_range),label=["bare series",nothing][2],alpha= 0.5,color = color_vec[1], order=1)

        #plot moment extrapolation (pade)
        x_range = collect(0:0.01:pade_on_moments_range_vec[m_idx])
        delta_ext_moment_pade = zeros(length(x_range))
        for (x_pos,x_val) in enumerate(x_range)
            δ_vec,r_vec = fromMomentsToδ([m(x_val) for m in m_vec_extrapolated_pade])
            delta_ext_moment_pade[x_pos] = δ_vec[m_idx]
        end
        plot!(plt_m_x[m_idx],x_range,delta_ext_moment_pade,label=["pade on moments",nothing][2],alpha= 0.5,linestyle=linestyle_vec[2],color = color_vec[2], order=1)

        #plot moment extrapolation (ida)
        if m_idx<5
            x_range = collect(0:0.01:ida_on_moments_range_vec[m_idx])
            delta_ext_moment_ida = zeros(length(x_range))
            for (x_pos,x_val) in enumerate(x_range)
                δ_vec,r_vec = fromMomentsToδ([m[round(Integer,1+1/0.01*x_val)] for m in m_vec_extrapolated_ida])
                delta_ext_moment_ida[x_pos] = δ_vec[m_idx]
            end
            plot!(plt_m_x[m_idx],x_range,delta_ext_moment_ida,label=["IDA on moments",nothing][2],alpha= 0.5,linestyle=linestyle_vec[3],color = color_vec[3], order=1)
        end

        #plot delta extrapolation (ida)
        if m_idx<5
            x_range = collect(0:0.01:ida_on_deltas_range_vec[m_idx])

            t = Taylor1(14-2*m_idx)

            taylor_exp = δ_vec_raw[m_idx](t)
            taylor_coefficients=[taylor_exp[i]  for i=0:length(taylor_exp)-1]

            for tay_idx =1:length(taylor_coefficients)
                if abs(taylor_coefficients[tay_idx]) < 0.0000000001
                    taylor_coefficients[tay_idx] = 0.0
                end
            end

            taylor_poly = Polynomial(taylor_coefficients)

            IDA_parameters=[2,5-m_idx,5-m_idx]

            IDA_approximant = get_intDiffApprox(taylor_poly,x_range,IDA_parameters[1],IDA_parameters[2],IDA_parameters[3])

            plot!(plt_m_x[m_idx],x_range,IDA_approximant,label=["IDA on deltas",nothing][2],alpha= 0.5,linestyle=linestyle_vec[4],color = color_vec[4], order=1)
        end

        #plot DMRG reference data 
        scatter!(plt_m_x[m_idx],betas,DMRG_deltas[k_pos][m_idx],label=["DMRG data",nothing][2],alpha= 0.5,color = color_vec[5], order=1)

    end
end

display(plot(plt_m_x...))








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


