##### Plotting script for the Pyrochlore lattice
# This script calculates and plots the dynamical structure factor for the S=1 Pyrochlore lattice using the DynHTE method and compares it to 
# INS data obtained and kindly provided by Plumb et al. Nature Physics volume 15, pages54–59 (2019)

###############################################
# Imports and Setup
###############################################

using JLD2,DelimitedFiles, HDF5, LaTeXStrings,Polynomials
using DynHTE

include("../plotConventions.jl")

###############################################
# Load DynHTE Data
###############################################
L = 12
spin_length = 1
hte_lattice = getLattice(L, "pyrochlore")

@load "examples/Pyrochlore_Lattice/Pyrochlore_Lattice_S1_c_iipDyn_nmax12_L12.jld2" c_iipDyn_mat
### c_iipDyn_mat can be loaded from a file, or calculated from scratch
### If you want to calculate it from scratch, uncomment the following lines:

#hte_graphs = load_dyn_hte_graphs(spin_length, L);
#c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice, hte_graphs);



###############################################
# Brillouin Zone Path Definitions
###############################################
#Avoiding the singularity at the Γ point
# Pyrochlore (2,2,l)
path1 = [(4pi, 4pi, -4pi + 0.01), (4pi, 4pi, 0), (4pi, 4pi, 4pi - 0.01)]
# Pyrochlore (h,h,2)
path2 = [(-4pi+ 0.01, -4pi+ 0.01, 4pi ), (0, 0, 4pi ), (4pi- 0.01, 4pi- 0.01, 4pi )]
pathticks = ["Γ", "W", "Γ"]

###############################################
# Calculation Parameters
###############################################
Nk = 56
n_max = L
r_max = 3
Jmev = 2.4
x0 = 4
w_vec = collect(0:0.055:12.01/Jmev)

###############################################
# Calculate Structure Factors
###############################################
function calculate_structure_factors(path)
    kvec, kticks_positions = create_brillouin_zone_path(path, Nk)
    m0_vec = Matrix{Float64}(undef, length(kvec), r_max+1)
    
    # Calculate moments
    for (kpos, k) in enumerate(kvec)
        c_kDyn = get_c_k(k, c_iipDyn_mat, hte_lattice)
        m_vec = get_moments_from_c_kDyn(c_kDyn)
        poly_x = Polynomial([0, 1], :x)
        
        for r in 0:r_max
            xm_norm_r = (m_vec[1+r] / m_vec[1+r](0)) * poly_x
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
    
    # Calculate JS matrix
    JS_mat = zeros(length(kvec), length(w_vec))
    for kpos in eachindex(kvec)
        δ_vec, r_vec = fromMomentsToδ(m0_vec[kpos,:])
        δ_vec_ext = extrapolate_δvec(δ_vec, length(δ_vec) - 1, length(δ_vec) - 1, 2000, true)
        JS_mat[kpos,:] = [JS(δ_vec_ext, 1.0 * x0, w, 0.02) for w in w_vec]
    end
    return JS_mat
end

###############################################
# Import and Process INS Data
###############################################
function importINSdata(filename)
    data = readdlm(filename)
    intensity = data[:, 1]
    err = data[:, 2]
    x = data[:, 3]
    y = data[:, 4]
    return intensity, err, x, y
end

# Process 22L data
filename = "examples/Pyrochlore_Lattice/E22l.txt"
intensity, err, x, y = importINSdata(filename)
x_unique_22l = sort(unique(x))
y_unique_22l = y[1:23]
z_22l = reshape(intensity, length(x_unique_22l), length(y_unique_22l))
err_22l = reshape(err, length(x_unique_22l), length(y_unique_22l))

# Process HH2 data
filename = "examples/Pyrochlore_Lattice/Ehh2.txt"
intensity, err, x, y = importINSdata(filename)
x_unique_hh2 = sort(unique(x))
indices = vcat([1,2], collect(5:19))
y_unique_hh2 = y[1:Int(length(intensity)/length(x_unique_hh2))][indices]
z_hh2 = reshape(intensity, length(x_unique_hh2), length(y[1:Int(length(intensity)/length(x_unique_hh2))]))[: ,indices]
err_hh2 = reshape(err, length(x_unique_hh2), length(y[1:Int(length(intensity)/length(x_unique_hh2))]))[: ,indices]

xINS = [x_unique_22l,x_unique_hh2]
yINS = [y_unique_22l,y_unique_hh2]
zINS = [z_22l,z_hh2]
errINS = [err_22l,err_hh2]

###############################################
# Plot Comparison Figures
###############################################
# Calculate structure factors
HTE_data = [calculate_structure_factors(path) for path in [path1, path2]]

# Plot functions
function plot_comparison_22l()

        title_string =["","","","","","","","",""] 
        
        xPlots,yPlots=2,2
        
        fig = CairoMakie.Figure(layout=(yPlots,xPlots), size=(aps_width,1.2*aps_width),fontsize=7)
        
        grid = fig[1,1]= GridLayout()
        
        plotsINS = []
        plotsDynHTE = []
        hm1 = 0
        hm2 = 0
        
        xlabels = [L"[22l]",L"[hh2]"]
        
        for i in 1:2
            
            ax=CairoMakie.Axis(grid[2,i],limits=(-2,2,0,12),xlabel=xlabels[i],ylabel=L"\omega [meV]",title=title_string[i],titlesize=10,xlabelsize=12,ylabelsize=12)
            hm1=CairoMakie.heatmap!(ax,xINS[i],yINS[i],zINS[i],colormap=:viridis,colorrange=(0.0,80),highclip=:white)
            push!(plotsINS, ax)
        
            #INS
            
          
            
            
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
        
        
        hideydecorations!(plotsINS[2])
        hidedecorations!(plotsDynHTE[2], grid = false)
        hidexdecorations!(plotsDynHTE[1])
        
        
        colgap!(subgrid, 0)
        rowgap!(subgrid, 0)
        
        colgap!(subgrid2, 0)
        rowgap!(subgrid2, 0)
        
        colgap!(grid, 7)
        rowgap!(grid, 8)
        
        
        resize_to_layout!(fig)
        display(fig)
        save("examples/Pyrochlore_Lattice/INS_comp.png",fig;px_per_unit=2.0)

        
        
end

function plot_full_comparison()
    #### PLot comparison of INS and DynHTE for 22L 
   
            
        title_string =["","","","","","","","",""] 
        
        xPlots,yPlots=2,1
        
        fig = CairoMakie.Figure(layout=(yPlots,xPlots), size=(aps_width,0.9*aps_width),fontsize=7, backgroundcolor = :transparent)
        
        grid = fig[1,1]= GridLayout()
        
        plotsINS = []
        plotsDynHTE = []
        plotsKslices = []
        hm1 = 0
        hm2 = 0
        
        xlabels = [L"[22l\,]",L"[hh2]"]
        
        for i in 1:1
            
            ax=CairoMakie.Axis(grid[2,2],limits=(-2,2,0,12),xlabel=xlabels[i],ylabel=L"\omega [meV]",title=title_string[i],titlesize=10,xlabelsize=12,ylabelsize=12)
            hm1=CairoMakie.heatmap!(ax,xINS[i],yINS[i],zINS[i],colormap=:viridis,colorrange=(0.0,80),highclip=:white)
            push!(plotsINS, ax)
        
            #INS
            
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
        cb2= CairoMakie.Colorbar(grid[1, 2],hm1,size=8,labelsize = 10,label = L"\frac{k_i \, \text{d}^2\sigma}{k_f \,\text{d}\Omega \text{d}E'} \;\;\; [\frac{\text{barn}}{\text{eV } \text{sr}}\text{ per Ni}]", vertical = false, )
        
        
         
        hideydecorations!(plotsINS[1], grid = false)
        
        
        colgap!(grid, 7)
        rowgap!(grid, 5)
        
        
        resize_to_layout!(fig)
        resize_to_layout!(fig)
    
        display(fig)
        save("examples/Pyrochlore_Lattice/INS_comp_22l.png",fig;px_per_unit = 600/96)
        

end

# Generate plots
plot_comparison_22l()
plot_full_comparison()
