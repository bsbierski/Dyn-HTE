using JLD2, LaTeXStrings
#using RobustPad

include("../../LatticeGraphs.jl")
include("../../Embedding.jl")
include("../../ConvenienceFunctions.jl")
include("../../plotConventions.jl")



L = 12

spin_length = 1/2
hte_graphs = load_dyn_hte_graphs(spin_length, L);
hte_lattice = getLattice(L, "pyrochlore");
@time c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice, hte_graphs);
@save "CaseStudy/Renormalized_Mean_Field_Theory/Pyrochlore_Lattice_S1half_c_iipDyn_nmax12_L12.jld2" c_iipDyn_mat


@load "CaseStudy/Pyrochlore_Lattice/Pyrochlore_Lattice_S1_c_iipDyn_nmax12_L12.jld2" c_iipDyn_mat
#2.Compute all correlations in the lattice
#@time Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);
@load "CaseStudy/Pyrochlore_Lattice/Correlation_DataS=S1half_L12.jld2" Correlators


#check the susceptibility with 10.1103/PhysRevB.53.14228
(brillouin_zone_cut([(0.0, 0.0, 0.0) (0.0, 0.0, 0.0); (0.0, 0.0, 0.0) (0.0, 0.0, 0.0)], Correlators, lattice, center_sites)[1]*4)[:, 1] .* [factorial(n) * 4^n for n in 0:max_order]

using Plots


### Compute A 2D Brillouin zone cut: 
N = 30
kx = range(-4pi, 4pi, length=N)
ky = range(-4pi, 4pi, length=N) #[0.] #for chains
kmat = [(x, x, y) for x in kx, y in ky]
c_kDyn =  get_c_k(kmat,c_iipDyn_mat,hte_lattice);

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
if true
    x = -2
    padetype = [6, 6]
    evaluate(y) = eval_correlator_LR_continuous_pad(y, x, padetype) #define evaluation function
    struc = (evaluate.(c_kDyn))
    p = Plots.heatmap(kx, ky, abs(x) * transpose(struc); clims=(0, 1), aspect_ratio=1 / sqrt(3), xlims=[-4pi, 4pi])
end






#Generate the path 

Nk = 62
kvec, kticks_positioins = create_brillouin_zone_path(path, Nk)



###### S(k,w) heatmap
using CairoMakie

x = 1.
w_vec = collect(0:0.314/2:6.0)
#JSkw_mat = get_JSkw_mat_finitex("total","padetanh",x,kvec,w_vec,0.02,3,3,200,false,Correlators,lattice,center_sites)
@time JSkw_mat = get_JSkw_mat("u_pade", x, kvec, w_vec, c_iipDyn_mat, hte_lattice; r_ext=200, f=0.4, r_min=3, r_max=3)
fig = Figure(fontsize=25, resolution=(400, 500));
ax = Axis(fig[1, 1], limits=(0, Nk + 1, 0, 6), ylabel=L"\omega/J=w", title="Pyrochlore S=1: x=" * string(x), titlesize=25, xlabelsize=25, ylabelsize=25, aspect=1 / 2);
hm = CairoMakie.heatmap!(ax, [k for k in 1:Nk+1], w_vec, JSkw_mat, colormap=:viridis, colorrange=(0.001, 0.3), highclip=:white);
ax.xticks = (kticks_positioins, pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm, size=40, label=L"J S(k,w)")
resize_to_layout!(fig);
display(fig)

using Plots
import Plots.plot!
###### dynamical Matsubara correlator (k-space)
n_max = 1 * L
if true
    k, k_label = (4pi, 4pi, 0*pi), "W"
    c_kDyn = get_c_k(k, c_iipDyn_mat, hte_lattice)
    m_vec = get_moments_from_c_kDyn(c_kDyn)
    poly_x = Polynomial([0, 1], :x)

    x_vec_bare = collect(0:0.025:1.4)
    x_vec = collect(0.0:0.1:4.0)


    ### with x-series

    plt_m = plot([0], [0], xlims=(0, x_vec[end]), ylims=(0, 0.5), label="", xlabel="x=J/T", ylabel=L"x \cdot m_r(x)/m_r(0)", legend=:topleft, title="Pyrochlore AFM S=1: moments @k=" * k_label * " x-pade")
    for r in 0:3
        #xm_norm_r = m_vec[1+r]
        xm_norm_r = (m_vec[1+r] / m_vec[1+r](0)) * poly_x
        println()
        println("r=$r")
        @show xm_norm_r
        plot!(plt_m, x_vec_bare, xm_norm_r.(x_vec_bare), color=color_vec[r+1], linewidth=0.4, label="r=$r")

        plot!(plt_m, x_vec, get_pade(xm_norm_r, 7 - r, 6 - r).(x_vec), color=color_vec[r+1], linestyle=linestyle_vec[2], label="")
        plot!(plt_m, x_vec, get_pade(xm_norm_r, 7 - r - 1, 6 - r - 1).(x_vec), color=color_vec[r+1], linestyle=linestyle_vec[3], label="")
        #   plot!(plt_m,x_vec,get_pade(xm_norm_r,6-r-2,6-r-2).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[4],label="")
    end
    xPlots, yPlots = 1, 1
    plt_final = plot(plt_m, layout=(yPlots, xPlots))
    #  display(plt_final)
    savefig(plt_final, "CaseStudy/Pyrochlore_Lattice/Pyrochlore_moments_x-series_k" * k_label * ".png")

    linestyle_vec[3]


    ### with u-series

    x0_vec = [0.5, 1.0, 2.0, 3.0]  # for these x the DSF will be computed
    r_max = 3         # maximal order for moment
    m0_vec = [Float64[] for _ in x0_vec]

    plt_m = plot([0], [0], xlims=(0, x_vec[end]), ylims=(0, 0.5), label="", xlabel=L"x=J/T", ylabel=L"x \cdot m_{\mathbf{k},2r}(x) \, / \, m_{\mathbf{k},2r}(0)", legend=:topleft, title="Pyrochlore AFM S=1: moments @k=" * k_label)
    plot!(plt_m, x_vec, -x_vec, color=:grey, label="x bare")
    plot!(plt_m, x_vec, -x_vec, color=:grey, linestyle=linestyle_vec[2], label="u Padé [7-r,6-r]")
    plot!(plt_m, x_vec, -x_vec, color=:grey, linestyle=linestyle_vec[3], label="u Padé [6-r,5-r]")
    annotate!(plt_m, 3, 2, Plots.text(L"\mathbf{k}=" * string(k_label), 7))

    #plot!(plt_m,title="Kagome AFM S=1/2: moment at k="*k_label*" (f=$f)")
    #plot!([2/f,2/f],[0,4],color=:grey,linewidth=10,label="u-freezing",alpha=0.3)
    for r in 0:r_max

        xm_norm_r = (m_vec[1+r] / m_vec[1+r](0)) * poly_x
        fun1 = get_pade(xm_norm_r, 7 - r, 6 - r)
        fun2 = get_pade(xm_norm_r, 6 - r, 5 - r)
        x_div = find_divergence_point(fun1, fun2, 0.002, x_min=0.0, x_max=x_vec[end], step=0.1)
        f = 1 / (2 * x_div)

        if  f < 0.4
            f = 0.3
        end
        
        if r == 0 && f >0.6
            f = 0.6
        end


        println(f)
        #  f= [0.3,0.3,0.5,1.9][r+1]
        u_vec = tanh.(f .* x_vec)
        u0_vec = tanh.(f .* x0_vec)


        xm_norm_r = coeffs(poly_x * m_vec[1+r] / (m_vec[1+r](0)))
        #= println()
        println("r=$r")
        @show xm_norm_r =#


        ufromx_mat = get_LinearTrafoToCoeffs_u(n_max + 1, f)
        p_u = Polynomial(ufromx_mat[1:n_max+2-2*r, 1:n_max+2-2*r] * xm_norm_r)

        plot!(plt_m, x_vec_bare, Polynomial(xm_norm_r).(x_vec_bare), color=color_vec[r+1], linewidth=0.4, label="r=$r", alpha=0.7)
        #plot!(plt_m,x_vec,p_u.(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[1],label="u-poly r=$r")

        plot!(plt_m, x_vec, get_pade(p_u, 7 - r, 6 - r).(u_vec), color=color_vec[r+1], linestyle=linestyle_vec[2], label="", alpha=0.7)
        plot!(plt_m, x_vec, get_pade(p_u, 6 - r, 5 - r).(u_vec), color=color_vec[r+1], linestyle=linestyle_vec[3], label="", alpha=0.7)
        #plot!(plt_m,x_vec,get_pade(p_u,5-r,4-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[4],label="",alpha=0.7)

        ### extract moments at x0_vec 
        for x0_pos in eachindex(x0_vec)
            x0 = x0_vec[x0_pos]
            u0 = u0_vec[x0_pos]
            push!(m0_vec[x0_pos], m_vec[1+r](0) / x0 * get_pade(p_u, 6 - r, 6 - r)(u0))
        end
    end

    xPlots, yPlots = 1, 1
    plt_final = plot(plt_m, layout=(yPlots, xPlots))
    #display(plt_final)
    savefig(plt_final, "CaseStudy/Pyrochlore_Lattice/Pyrochlore_moments_u-series_k" * k_label * ".png")




    ### DSF and δ_r (as inset) for x ∈ x0_vec
    w_vec = collect(0.0:0.025:4)
    plt_JS = plot(color=:grey, legend=:bottomleft, label="Dyn-HTE", xlims=(0, w_vec[end]), ylims=(0, 0.6), xlabel=L"\omega/J=w", ylabel=L"J\, S(\mathbf{k}=" * k_label * L",\omega)")

    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]

        ### plot Dyn-HTE
        δ_vec, r_vec = fromMomentsToδ(m0_vec[x0_pos])
        #scatter!(plt_δ,r_vec,δ_vec,color=thermalCol4_vec[x0_pos],label="x=$x0")
        δ_vec_ext = extrapolate_δvec(δ_vec, length(δ_vec) - 1, length(δ_vec) - 1, 1000, true)
        #plot!(plt_δ,3:7,δ_vec_ext[4:8],label="",color=thermalCol4_vec[x0_pos]) ###test extrapolation of δr
        plot!(plt_JS, w_vec, [JS(δ_vec_ext, 1.0 * x0, w, 0.02) for w in w_vec], color=thermalCol4_vec[x0_pos], label="")
    end
    xPlots, yPlots = 1, 1
    plt_final = plot(plt_JS)
    #display(plt_final)
    savefig(plt_final, "CaseStudy/Pyrochlore_Lattice/Pyrochlore_JS_k" * k_label * ".png")

    println("done")
end


xINS[1]



############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone
#---pyrochlore (2,2, l)
path1 = [(4pi, 4pi, -4pi + 0.01), (4pi, 4pi, 0), (4pi, 4pi, 4pi - 0.01)]
pathticks = ["Γ", "W", "Γ"]

#---pyrochlore (h,h, 2)
path2 = [(-4pi+ 0.01, -4pi+ 0.01, 4pi ), (0, 0, 4pi ), (4pi- 0.01, 4pi- 0.01, 4pi )]
pathticks = ["Γ", "W", "Γ"]
HTE_data = []

Nk = 56
n_max = L
kvec, kticks_positioins = create_brillouin_zone_path(path1, Nk)
r_max = 3
m0_vec = Matrix{Float64}(undef,length(kvec),r_max+1)
x0 = 1
w_vec = []

if true 
HTE_data = []
Jmev = 2.4
for path in [path1,path2]


Nk = 56
n_max = L
kvec, kticks_positioins = create_brillouin_zone_path(path, Nk)
r_max = 3
m0_vec = Matrix{Float64}(undef,length(kvec),r_max+1)
x0 = 1



for (kpos, k) in enumerate(kvec)

    c_kDyn = get_c_k(k, c_iipDyn_mat, hte_lattice)
    m_vec = get_moments_from_c_kDyn(c_kDyn)
    poly_x = Polynomial([0, 1], :x)

    for r in 0:r_max

        xm_norm_r = (m_vec[1+r] / m_vec[1+r](0)) * poly_x
        fun1 = get_pade(xm_norm_r, 7 - r, 6 - r)
        fun2 = get_pade(xm_norm_r, 6 - r, 5 - r)
        x_div = find_divergence_point(fun1, fun2, 0.002, x_min=0.0, x_max=4, step=0.1)
        f = 1 / (2 * x_div)

        
        if r == 0 && f >0.6
            f = 0.6
        end

        f = 0.6

        if f > 0.3
            xm_norm_r = coeffs(poly_x * m_vec[1+r] / (m_vec[1+r](0)))

            ufromx_mat = get_LinearTrafoToCoeffs_u(n_max + 1, f)
            p_u = Polynomial(ufromx_mat[1:n_max+2-2*r, 1:n_max+2-2*r] * xm_norm_r)
    
            u0 = tanh.(f .* x0)
            m0_vec[kpos, r+1] = m_vec[1+r](0)/x0 * get_pade(p_u, 6 - r, 6 - r)(u0)
        
        else
            xm_norm_r = coeffs(poly_x * m_vec[1+r] / (m_vec[1+r](0)))
            p_x = Polynomial(xm_norm_r)
            m0_vec[kpos, r+1] = m_vec[1+r](0)/x0 * get_pade(p_x, 6 - r, 6 - r)(x0)
        end

    end

end


w_vec = collect(-12/Jmev:0.055:12.01/Jmev)
JS_mat = zeros(length(kvec),length(w_vec))
for kpos in eachindex(kvec)
        δ_vec, r_vec = fromMomentsToδ(m0_vec[kpos,:])
        #scatter!(plt_δ,r_vec,δ_vec,color=thermalCol4_vec[x0_pos],label="x=$x0")
        δ_vec_ext = extrapolate_δvec(δ_vec, length(δ_vec) - 1, length(δ_vec) - 1, 200, true)
        JS_mat[kpos,:] = [#= exp(-0.1*kvec[kpos]*kvec[kpos])* =#JS(δ_vec_ext, 1.0 * x0, w, 0.02) for w in w_vec]
end

plt = Plots.heatmap([1:length(kvec)],w_vec*Jmev,transpose(JS_mat); clims=(0, 0.3), size = (260,400) , colorbar = true, colormap = :viridis)
savefig(plt, "CaseStudy/Pyrochlore_Lattice/test.png")



push!(HTE_data, JS_mat)

end
end


if true

# Load necessary package
using DelimitedFiles


# Read the data from the file



# Extract columns into separate vectors
function importINSdata(filename)
data = readdlm(filename)
intensity = data[:, 1]
err     = data[:, 2]
x         = data[:, 3]
y         = data[:, 4]

return (intensity,err,x,y)
end
    xINS = []
    yINS = []
    zINS = []
    errINS = []
    

####22L plot
# Define the path to your file
filename = "CaseStudy/Pyrochlore_Lattice/E22l.txt"  # Replace with your actual file path
(intensity,err,x,y) = importINSdata(filename)
# Create a grid from x and y values (assumes regular grid, sorted)
x_unique = sort(unique(x))
y_unique = y[1:23]
# Reshape intensity to match the grid dimensions
z = reshape(intensity,  length(x_unique),length(y_unique))
err = reshape(err,  length(x_unique),length(y_unique))
# Plot heatmap
Plots.heatmap(x_unique,y_unique,z',xlabel="[22L]", ylabel="ω(meV)", colorbar_title="Intensity", title="INS-data", clims=(0,80.0), xlims = (-2,2.),colormap = :viridis)

push!(xINS, x_unique)
push!(yINS, y_unique)
push!(zINS, z)
push!(errINS, err)

filename = "CaseStudy/Pyrochlore_Lattice/Ehh2.txt"  # Replace with your actual file path
(intensity,err,x,y) = importINSdata(filename)
x_unique = sort(unique(x))
y_unique = y[1:Int(length(intensity)/length(x_unique))]

z = reshape(intensity,  length(x_unique),length(y_unique))
err = reshape(err,  length(x_unique),length(y_unique))
# Plot heatmap
indices = vcat([1,2],collect(5:19))

Plots.heatmap(x_unique,y_unique,z',xlabel="[HH2]", ylabel="ω(meV)", colorbar_title="Intensity", title="INS-data", clims=(0,80.0), xlims = (-2,2.),colormap = :viridis)

push!(xINS, x_unique)
push!(yINS, y_unique[indices])
push!(zINS, z[:,indices])
push!(errINS, err[:,indices])

end 



####PLotting
if true
kvec

betas = [0,2]

title_string =["","","","","","","","",""] 

xPlots,yPlots=2,2

fig = CairoMakie.Figure(layout=(yPlots,xPlots), size=(aps_width,1.2*aps_width),fontsize=7)

grid = fig[1,1]= GridLayout()

plotsDMRG = []
plotsDynHTE = []
plotsKslices = []
hm1 = 0
hm2 = 0

xlabels = [L"[22l]",L"[hh2]"]

for (i,beta) in enumerate(betas)
    
    ax=CairoMakie.Axis(grid[2,i],limits=(-2,2,0,12),xlabel=xlabels[i],ylabel=L"\omega [meV]",title=title_string[i],titlesize=10,xlabelsize=12,ylabelsize=12)
    hm1=CairoMakie.heatmap!(ax,xINS[i],yINS[i],zINS[i],colormap=:viridis,colorrange=(0.0,80),highclip=:white)
    push!(plotsDMRG, ax)

    #DMRG
    
  
    
    
    #dynHTE spot
    ax_dynHTE=CairoMakie.Axis(grid[1,i],limits=(-2,2,0,12),ylabel=L"\omega [meV]",title=title_string[i+3],titlesize=10,xlabelsize=12,ylabelsize=12)
    hm2=CairoMakie.heatmap!(ax_dynHTE,-2:4/Nk:2,w_vec*Jmev,HTE_data[i],colormap=:viridis,colorrange=(0.0,0.45),highclip=:white)
    push!(plotsDynHTE, ax_dynHTE)

    ###annotate
    CairoMakie.text!(ax_dynHTE,  -1.8, 11.7, text=["Dyn-HTE",""][i], fontsize=10, color=:white, align=(:left, :top) )
    CairoMakie.text!(ax_dynHTE, -1.8, 11.7, text=["","x = $(Int(x0))"][i], fontsize=10, color=:white, align=(:left, :top) )
    CairoMakie.text!(ax,  -1.8, 11.7, text=["INS",""][i], fontsize=10, color=:white, align=(:left, :top) )

    CairoMakie.text!(ax, 1.3, 11.7, text=["c)","d)"][i], fontsize=10, color=:white, align=(:left, :top) )
    CairoMakie.text!(ax_dynHTE, 1.3, 11.7, text=["a)","b)"][i], fontsize=10, color=:white, align=(:left, :top) )

end 

subgrid = GridLayout(grid[1, 3], tellheight = false)
subgrid2 = GridLayout(grid[2, 3], tellheight = false)

Label(subgrid[1, 1], L"JS(k,\omega)",fontsize=10)
Label(subgrid2[1, 1],L"\frac{k_i \text{d}^2\sigma}{k_f \text{d}\Omega \text{d}E'}",fontsize=10)


cb1=CairoMakie.Colorbar(subgrid[2, 1],hm2,size=11,labelsize = 10 , colorrange = (0,1), ) 
cb2= CairoMakie.Colorbar(subgrid2[2, 1],hm1,size=11,labelsize = 10,label = L"\text{barn eV}^{-1} \text{sr}^{-1} \text{ per Ni}")


hideydecorations!(plotsDMRG[2])
hidedecorations!(plotsDynHTE[2], grid = false)
hidexdecorations!(plotsDynHTE[1])


colgap!(subgrid, 0)
rowgap!(subgrid, 0)

colgap!(subgrid2, 0)
rowgap!(subgrid2, 0)

colgap!(grid, 7)
rowgap!(grid, 8)


resize_to_layout!(fig)
# Plots.savefig("dynHTE_dmrg_test.pdf")
display(fig)
save("CaseStudy/Pyrochlore_Lattice/INS_comp.png",fig;px_per_unit=2.0)

end

if true
    kvec
    
    betas = [0]
    
    title_string =["","","","","","","","",""] 
    
    xPlots,yPlots=2,1
    
    fig = CairoMakie.Figure(layout=(yPlots,xPlots), size=(aps_width,0.9*aps_width),fontsize=7, backgroundcolor = :transparent)
    
    grid = fig[1,1]= GridLayout()
    
    plotsDMRG = []
    plotsDynHTE = []
    plotsKslices = []
    hm1 = 0
    hm2 = 0
    
    xlabels = [L"[22l\,]",L"[hh2]"]
    
    for (i,beta) in enumerate(betas)
        
        ax=CairoMakie.Axis(grid[2,2],limits=(-2,2,0,12),xlabel=xlabels[i],ylabel=L"\omega [meV]",title=title_string[i],titlesize=10,xlabelsize=12,ylabelsize=12)
        hm1=CairoMakie.heatmap!(ax,xINS[i],yINS[i],zINS[i],colormap=:viridis,colorrange=(0.0,80),highclip=:white)
        push!(plotsDMRG, ax)
    
        #DMRG
        
        #dynHTE spot
        ax_dynHTE=CairoMakie.Axis(grid[2,1],limits=(-2,2,0,12),xlabel=xlabels[i],ylabel=L"\omega [meV]",title=title_string[i+3],titlesize=10,xlabelsize=12,ylabelsize=12)
        hm2=CairoMakie.heatmap!(ax_dynHTE,-2:4/Nk:2,w_vec*Jmev,HTE_data[i],colormap=:viridis,colorrange=(0.0,0.45),highclip=:white)
        push!(plotsDynHTE, ax_dynHTE)
    
        ###annotate
        CairoMakie.text!(ax_dynHTE,  -1.8, 11.7, text=["Dyn-HTE",""][i], fontsize=10, color=:white, align=(:left, :top) )
        CairoMakie.text!(ax_dynHTE, -1.6, 10.9, text=["x = $(Int(x0))"][i], fontsize=10, color=:white, align=(:left, :top) )
        CairoMakie.text!(ax,  -1.8, 11.7, text=["INS",""][i], fontsize=10, color=:white, align=(:left, :top) )
    
        CairoMakie.text!(ax, 1.6, 11.7, text=["b)"][i], fontsize=10, color=:white, align=(:left, :top) )
        CairoMakie.text!(ax_dynHTE, 1.6, 11.7, text=["a)"][i], fontsize=10, color=:white, align=(:left, :top) )
    
    end 
    
    cb1=CairoMakie.Colorbar(grid[1, 1],hm2,size=10,labelsize = 10 , label = L"JS(\mathbf{k},\omega)", colorrange = (0,1), vertical = false, ) 
    cb2= CairoMakie.Colorbar(grid[1, 2],hm2,size=8,labelsize = 10,label = L"\frac{k_i \, \text{d}^2\sigma}{k_f \,\text{d}\Omega \text{d}E'} \;\;\; [\frac{\text{barn}}{\text{eV } \text{sr}}\text{ per Ni}]", vertical = false, )
    
    
     
    hideydecorations!(plotsDMRG[1], grid = false)
    
    
    colgap!(subgrid, 0)
    rowgap!(subgrid, 0)
    
    colgap!(subgrid2, 0)
    rowgap!(subgrid2, 0)
    
    colgap!(grid, 7)
    rowgap!(grid, 5)
    
    
    resize_to_layout!(fig)
    resize_to_layout!(fig)

    # Plots.savefig("dynHTE_dmrg_test.pdf")
    display(fig)
    save("CaseStudy/Pyrochlore_Lattice/INS_comp_22l.png",fig;px_per_unit = 600/96)
    
end

aps_width




if true 
    scale = 360
    Jmev = 2.4
plot(yINS[1],zINS[1][61,:]/scale);
plot!(yINS[1],zINS[1][60,:]/scale);
plot!(w_vec*Jmev,HTE_data[1][29,:],);


plot!(yINS[1],zINS[1][65,:]/scale);
plot!(yINS[1],zINS[1][66,:]/scale);
plot!(w_vec*Jmev,HTE_data[1][36,:],)

plot!(yINS[1],zINS[1][70,:]/scale);
plot!(yINS[1],zINS[1][71,:]/scale);
plot!(w_vec*Jmev,HTE_data[1][43,:],)



end


######## CONSTANT ENERGY


if true 
    HTE_data = []
    Jmev = 2.4
    for path in [path1,path2]
    
    
    Nk = 56
    n_max = L
    kvec, kticks_positioins = create_brillouin_zone_path(path, Nk)
    r_max = 3
    m0_vec = Matrix{Float64}(undef,length(kvec),r_max+1)
    x0 = 0.4
    
    
    
    for (kpos, k) in enumerate(kvec)
    
        c_kDyn = get_c_k(k, c_iipDyn_mat, hte_lattice)
        m_vec = get_moments_from_c_kDyn(c_kDyn)
        poly_x = Polynomial([0, 1], :x)
    
        for r in 0:r_max
    
            xm_norm_r = (m_vec[1+r] / m_vec[1+r](0)) * poly_x
            fun1 = get_pade(xm_norm_r, 7 - r, 6 - r)
            fun2 = get_pade(xm_norm_r, 6 - r, 5 - r)
            x_div = find_divergence_point(fun1, fun2, 0.002, x_min=0.0, x_max=4, step=0.1)
            f = 1 / (2 * x_div)
    
            
            if r == 0 && f >0.6
                f = 0.6
            end
    
            f = 0.6
    
            if f > 0.3
                xm_norm_r = coeffs(poly_x * m_vec[1+r] / (m_vec[1+r](0)))
    
                ufromx_mat = get_LinearTrafoToCoeffs_u(n_max + 1, f)
                p_u = Polynomial(ufromx_mat[1:n_max+2-2*r, 1:n_max+2-2*r] * xm_norm_r)
        
                u0 = tanh.(f .* x0)
                m0_vec[kpos, r+1] = m_vec[1+r](0)/x0 * get_pade(p_u, 6 - r, 6 - r)(u0)
            
            else
                xm_norm_r = coeffs(poly_x * m_vec[1+r] / (m_vec[1+r](0)))
                p_x = Polynomial(xm_norm_r)
                m0_vec[kpos, r+1] = m_vec[1+r](0)/x0 * get_pade(p_x, 6 - r, 6 - r)(x0)
            end
    
        end
    
    end
    
    
    w_vec = unique(yINS[1])/Jmev
    JS_mat = zeros(length(kvec),length(w_vec))
    for kpos in eachindex(kvec)
            δ_vec, r_vec = fromMomentsToδ(m0_vec[kpos,:])
            #scatter!(plt_δ,r_vec,δ_vec,color=thermalCol4_vec[x0_pos],label="x=$x0")
            δ_vec_ext = extrapolate_δvec(δ_vec, length(δ_vec) - 1, length(δ_vec) - 1, 600, true)
            JS_mat[kpos,:] = [JS(δ_vec_ext, 1.0 * x0, w, 0.02) for w in w_vec]
    end
    
    plt = heatmap([1:length(kvec)],w_vec*Jmev,transpose(JS_mat); clims=(0, 0.3), size = (260,400) , colorbar = true, colormap = :viridis)
    savefig(plt, "CaseStudy/Pyrochlore_Lattice/test.png")
    
    push!(HTE_data, JS_mat)
    
    end
    end




if true 
    scale = 420

    plt = plot(-2:4/Nk:2,HTE_data[1][:,6], xlims = (-2,2), ylims = (0,0.2), color=color_vec[1], label=L"\omega="* string(yINS[1][6]), size=(aps_width*xPlots,0.62*aps_width*yPlots))
    scatter!(xINS[1],zINS[1][:,6]/scale, color=color_vec[1], label = "")
    plot!(-2:4/Nk:2,HTE_data[1][:,10], color=color_vec[2], label=L"\omega="* string(yINS[1][10]))
    scatter!(xINS[1],zINS[1][:,10]/scale, yerror  = errINS[1][:,10]/scale, color=color_vec[2], label = "")
    plot!(-2:4/Nk:2,HTE_data[1][:,13], color=color_vec[3], label=L"\omega="* string(yINS[1][13]))
    scatter!(xINS[1],zINS[1][:,13]/scale, yerror  = errINS[1][:,13]/scale, color=color_vec[3], label = "")
    plot!(-2:4/Nk:2,HTE_data[1][:,17], color=color_vec[4], label=L"\omega="* string(yINS[1][17]))
    scatter!(xINS[1],zINS[1][:,17]/scale , yerror  = errINS[1][:,17]/scale , color=color_vec[4], label = "")
    plot!([],[],color = :black, label = "Dyn-HTE")
    scatter!([],[],color = :black, label = "INS")
    

    savefig(plt, "CaseStudy/Pyrochlore_Lattice/EqualFrequency.png")
end

using Plots
gr()

x = 1:5
y = rand(5)
yerr = [0.1, 0.2, 0.15, 0.05, 0.3]

plot(x, y, err = yerr, seriestype = :scatter, label = "")



# Scatter plot with error bars
 plot(x, y; yerror = yerr, label = "Scatter with Error Bars", marker = :circle, legend = :top)



######## CONSTANT k



if true 
    HTE_data = []
    Jmev = 2.8
    for path in [path1,path2]
    
    
    Nk = 56
    n_max = L
    kvec, kticks_positioins = create_brillouin_zone_path(path, Nk)
    r_max = 3
    m0_vec = Matrix{Float64}(undef,length(kvec),r_max+1)
    x0 = 0.5
    
    
    
    for (kpos, k) in enumerate(kvec)
    
        c_kDyn = get_c_k(k, c_iipDyn_mat, hte_lattice)
        m_vec = get_moments_from_c_kDyn(c_kDyn)
        poly_x = Polynomial([0, 1], :x)
    
        for r in 0:r_max
    
            xm_norm_r = (m_vec[1+r] / m_vec[1+r](0)) * poly_x
            fun1 = get_pade(xm_norm_r, 7 - r, 6 - r)
            fun2 = get_pade(xm_norm_r, 6 - r, 5 - r)
            x_div = find_divergence_point(fun1, fun2, 0.002, x_min=0.0, x_max=4, step=0.1)
            f = 1 / (2 * x_div)
    
            
            if r == 0 && f >0.6
                f = 0.6
            end
    
            f = 0.6
    
            if f > 0.3
                xm_norm_r = coeffs(poly_x * m_vec[1+r] / (m_vec[1+r](0)))
    
                ufromx_mat = get_LinearTrafoToCoeffs_u(n_max + 1, f)
                p_u = Polynomial(ufromx_mat[1:n_max+2-2*r, 1:n_max+2-2*r] * xm_norm_r)
        
                u0 = tanh.(f .* x0)
                m0_vec[kpos, r+1] = m_vec[1+r](0)/x0 * get_pade(p_u, 6 - r, 6 - r)(u0)
            
            else
                xm_norm_r = coeffs(poly_x * m_vec[1+r] / (m_vec[1+r](0)))
                p_x = Polynomial(xm_norm_r)
                m0_vec[kpos, r+1] = m_vec[1+r](0)/x0 * get_pade(p_x, 6 - r, 6 - r)(x0)
            end
    
        end
    
    end
    
    
    w_vec = unique(yINS[1])/Jmev
    JS_mat = zeros(length(kvec),length(w_vec))
    for kpos in eachindex(kvec)
            δ_vec, r_vec = fromMomentsToδ(m0_vec[kpos,:])
            #scatter!(plt_δ,r_vec,δ_vec,color=thermalCol4_vec[x0_pos],label="x=$x0")
            δ_vec_ext = extrapolate_δvec(δ_vec, length(δ_vec) - 1, length(δ_vec) - 1, 600, true)
            JS_mat[kpos,:] = [JS(δ_vec_ext, 1.0 * x0, w, 0.02) for w in w_vec]
    end
    
    plt = heatmap([1:length(kvec)],w_vec*Jmev,transpose(JS_mat); clims=(0, 0.3), size = (260,400) , colorbar = true, colormap = :viridis)
    savefig(plt, "CaseStudy/Pyrochlore_Lattice/test.png")
    

    push!(HTE_data, JS_mat)
    
    end
    end



if true 
    scale = 370
    Jmev = 2.4
plt = scatter(yINS[1],zINS[1][61,:]/scale,  yerr = errINS[1][61,:]/scale , xlims = (0.4,12),  ylims = (0.0,0.6),color=color_vec[1] , label  = "hh0", title = "x = " *string(x0));
scatter!(yINS[1],zINS[1][60,:]/scale,color=color_vec[1], label = "");
plot!(w_vec*Jmev,HTE_data[1][29,:],color=color_vec[1], label = "");


scatter!(yINS[1],zINS[1][65,:]/scale,color=color_vec[2], label  = "hh1/2");
scatter!(yINS[1],zINS[1][66,:]/scale,color=color_vec[2], label = "");
plot!(w_vec*Jmev,HTE_data[1][36,:],color=color_vec[2], label = "")

scatter!(yINS[1],zINS[1][70,:]/scale,color=color_vec[3], label  = "hh1");
scatter!(yINS[1],zINS[1][71,:]/scale,color=color_vec[3], label = "");
plot!(w_vec*Jmev,HTE_data[1][43,:],color=color_vec[3], label = "")

savefig(plt, "CaseStudy/Pyrochlore_Lattice/Equalk.png")


end
plt