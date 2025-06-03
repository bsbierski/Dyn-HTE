###############################################
# Imports and Setup
###############################################
using DynHTE
using CairoMakie
using Symbolics, TaylorSeries,Polynomials

###############################################
# Load Lattice Data
###############################################
L = 12
chain_lattice = getLattice(L,"chain");
triang_lattice = getLattice(L,"triang");
kagome_lattice = getLattice(L,"kagome");
pyrochlore_lattice = getLattice(L,"pyrochlore");

###############################################
# Load Correlation Data
###############################################
@load "examples/Renormalized_Mean_Field_Theory/Chain_S1half_c_iipDyn_nmax12_L12.jld2" c_iipDyn_mat
c_iip_chain = c_iipDyn_mat

@load "examples/Renormalized_Mean_Field_Theory/Triangular_Lattice_S1half_c_iipDyn_nmax12_L12.jld2" single_stored_object
c_iip_triang = single_stored_object

@load "examples/Renormalized_Mean_Field_Theory/Kagome_Lattice_S1half_c_iipDyn_nmax12_L12.jld2" c_iipDyn_mat
c_iip_kagome = c_iipDyn_mat

@load "examples/Renormalized_Mean_Field_Theory/Pyrochlore_Lattice_S1half_c_iipDyn_nmax12_L12.jld2" c_iipDyn_mat
c_iip_pyrochlore = c_iipDyn_mat

###############################################
# QMC Reference Data
###############################################
# QMC data generated with the QMC Worm Code: N. Sadoune and L. Pollet, https://github.com/lodepollet/worm (2024).
QMCbetas = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 4, 5, 6, 7, 8, 9, 10]
fs = [4, 4.228359252922975, 4.972941835020016, 6.262842981660236, 
  8.034188800241, 10.170915326395479, 12.513608854443149, 
  17.500654887074752, 22.59466799494428, 27.647081475999613, 
  32.696822363289336, 37.66727122510188, 42.631499009231305, 
  47.63892189825918]
gs = [0, 0.5331037777143603, 1.1749560328465218, 1.9830297440003888, 
  2.9457772366259505, 4.039309023432947, 5.17629294244374, 
  7.486534321669895, 9.748107036203425, 11.935348427190947, 
  14.110065180402612, 16.219086068396088, 18.32277163971843, 
  20.458603889181074]
ϵ2s = [0, 0.0003315853417177263, 0.0019288270029839015, 
  0.00799150745186226, 0.031118623522432845, 0.06847416640854603, 
  0.13502108782626027, 0.3125184265080508, 0.5122626846871366, 
  0.7048405267364483, 0.8676720271255824, 1.030294622644043, 
  1.1755097262206777, 1.3017551127324813]
ϵ3s = [0, 0.0003993867686287786, 0.0020677244995818476, 0.00048139162916391624, 
  0.000951201316996593, 0.004286699078212941, 0.01780962619837747, 
  0.0600786797894247, 0.14557215813130953, 0.25765708894019, 
  0.37717518289526825, 0.49576783668838303, 0.6078040978672965, 
  0.7371619522846239]
ϵ4s = [0, 0.001503856603280766, 0.0017316957332395109, 0.0031829822967957945, 
  0.005804781686516066, 0.0002142847371912503, 2.9629971502247966e-05, 
  0.01477502227527624, 0.03337474004648612, 0.06810672839164181, 
  0.1104637284387855, 0.15535103083658305, 0.21985050785598773, 
  0.25459888419265553]

###############################################
# Helper Functions
###############################################
"""
    get_unique_distances(lattice)
Returns tuple of unique distances and corresponding site indices from reference site in lattice.
"""
function get_unique_distances(lattice)
    site_pos = lattice.lattice.sitePositions
    ref_site_pos = lattice.lattice.sitePositions[lattice.basis_positions[1]]
    
    # Calculate all distances from reference site
    distances = sqrt.(dot.(site_pos .- Ref(ref_site_pos), site_pos .- Ref(ref_site_pos)))
    
    # Find unique distances with numerical tolerance
    unique_distances = Float64[]
    indices_by_distance = Vector{Int}[]
    for (i, d) in enumerate(distances)
        found = false
        for (j, ud) in enumerate(unique_distances)
            if isapprox(d, ud, rtol=1e-10)
                push!(indices_by_distance[j], i)
                found = true
                break
            end
        end
        if !found
            push!(unique_distances, d)
            push!(indices_by_distance, [i])
        end
    end
    
    # Sort by distance
    p = sortperm(unique_distances)
    
    return unique_distances[p], indices_by_distance[p]
end

function calc_taylorinvmat_fun(corr)
    orderJ,orderω = size(corr)
    @variables x
    taylormat = y -> sum(y[i,1]*x^(i-1) for i =  1:orderJ)
    invtaylormat = inv(taylormat(corr));
    t = Taylor1(Float64,orderJ-1)
    return  substitute(invtaylormat, Dict(x=>t)) 
end

import Base
function Base.Float64(x::Num)
    return Float64(Symbolics.value(x))
end

function calc_taylorinvmat_fun(corr::Matrix{Matrix{Float64}})::Matrix{Taylor1{Float64}}
    orderJ,orderω = size(corr[1])
    @variables x::Real
    taylormat = y -> sum(y[i,1]*x^(i-1) for i =  1:orderJ)
    invtaylormat = inv(taylormat.(corr));
    t = Taylor1(Float64,orderJ-1)
    subsfun = y -> substitute(real(y), Dict(x=>t)) 
    taylorinvmat_fun =  subsfun.(invtaylormat)
    return taylorinvmat_fun
end

###############################################
# Calculate Correlations
###############################################
# Chain lattice calculations
N = 32
kx = (1:N)*2pi/(N)
kmat = Tuple.([[x,y] for x in kx , y in [0]])
c_kDyn = get_c_k(kmat,c_iip_chain,chain_lattice)
invstruc = calc_taylorinvmat_fun.(c_kDyn)
invcorrs_chain = inverse_fourier_transform(kmat,invstruc,chain_lattice)

# Triangular lattice calculations
N = 50
kx = (1:N)*2pi/(N)
ky = (1:N)*2pi/(N)
kmat = Tuple.([[1,1/sqrt(3)].*x .+ [1,-1/sqrt(3)].*y for x in kx, y in ky])
c_kDyn = get_c_k(kmat,c_iip_triang,triang_lattice)
invstruc = calc_taylorinvmat_fun.(c_kDyn)
invcorrs_triang = inverse_fourier_transform(kmat,invstruc,triang_lattice)

# Kagome lattice calculations
N = 30
kx = (1:N)*2pi/(N)
ky = (1:N)*2pi/(N)
kmat = Tuple.([[1,1/sqrt(3)].*x .+ [1,-1/sqrt(3)].*y for x in kx, y in ky])
c_kDyn = get_c_k_subl(kmat,c_iip_kagome,kagome_lattice)
invstruc = Array{Matrix{Taylor1{Float64}}}(undef,N,N)
for i=1:N, j=1:N
    invstruc[i,j] = calc_taylorinvmat_fun(c_kDyn[i,j])
end
invcorrs_kagome = inverse_fourier_transform_subl(kmat,invstruc,kagome_lattice);


###### Pyrochlore lattice calculations: Take very long.

#= 
N = 26;
kx = (1:N)*2pi/(N);
ky = (1:N)*2pi/(N);
kz = (1:N)*2pi/(N); 
kmat = Tuple.([[-1,1,1].*x .+ [1,-1,1].*y .+ [1,1,-1].*z for x in kx, y in ky , z in kz ]);
@time c_kDyn =  get_c_k_subl(kmat,c_iip_pyrochlore,pyrochlore_lattice);
#test if the inverse is correct
invcorrstest = inverse_fourier_transform_subl(kmat,c_kDyn,pyrochlore_lattice);
maximum(maximum.(invcorrstest.-c_iip_pyrochlore))

#calculate f and g
invstruc = Array{Matrix{Taylor1{Float64}}}(undef,N,N,N)
#= Threads.@threads =# for i=1:N
    for j = 1:N
      for k = 1:N
        invstruc[i,j,k] = calc_taylorinvmat_fun(c_kDyn[i,j,k])

    end
end
end
invcorrs_pyrochlore = inverse_fourier_transform_subl(kmat,invstruc,pyrochlore_lattice); =#
@load "examples/Renormalized_Mean_Field_Theory/rMFparams_pyrochlore.jld2" invcorrs_pyrochlore

###############################################
# Plotting
###############################################
if true
    # Create figure and axes
    fig = CairoMakie.Figure(layout=(3,3), size=(aps_width,aps_width),fontsize=7)
    grid = fig[1,1]= GridLayout()
    axchain = CairoMakie.Axis(grid[1,1], title = "Chain" ,xlabel=L"x = J/T" ,titlesize=10,xlabelsize=12,ylabelsize=12,xgridvisible = false , limits=(0,6,-0.05,1.05), yticks = [0,0.2,0.4,0.6,0.8,1.0])
    axtriang = CairoMakie.Axis(grid[1,2], title = "Triangular Lattice",xlabel=L"x = J/T" ,titlesize=10,xlabelsize=12,ylabelsize=12,xgridvisible = false, limits=(0,6,-0.05,1.05), yticks = [0,0.2,0.4,0.6,0.8,1.0])
    axkagome = CairoMakie.Axis(grid[2,1], title = "Kagome Lattice", xlabel=L"x = J/T" ,titlesize=10,xlabelsize=12,ylabelsize=12,xgridvisible = false, limits=(0,6,-0.05,1.05), yticks = [0,0.2,0.4,0.6,0.8,1.0])
    axpyro = CairoMakie.Axis(grid[2,2], title = "Pyrochlore Lattice", xlabel=L"x = J/T" ,titlesize=10,xlabelsize=12,ylabelsize=12,xgridvisible = false, limits=(0,6,-0.05,1.05), yticks = [0,0.2,0.4,0.6,0.8,1.0])

    # Plot settings
    xvec = collect(0:0.1:10)

    # Chain lattice plot
    if true
        f1 = 0.4
        f2 = 0.4
        f3 = 0.4
        transf1 = get_LinearTrafoToCoeffs_u(12,f1)
        transf2 = get_LinearTrafoToCoeffs_u(13,f2)
        transf3 = get_LinearTrafoToCoeffs_u(12,f3)
        dists, indices = get_unique_distances(chain_lattice)
        invcorrs = invcorrs_chain
        gbyf = transf1*(invcorrs[indices[2][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ2byf = transf3*(invcorrs[indices[3][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ3byf = transf3*(invcorrs[indices[4][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ4byf = transf3*(invcorrs[indices[5][1]]/invcorrs[indices[1][1]]).coeffs
        byf = transf2*vcat(0,(1/invcorrs[indices[1][1]]).coeffs)[:]
        
        pade_gbyf = get_pade(Polynomial(gbyf),6,6)
        pade_ϵ2byf = get_pade(Polynomial(ϵ2byf),6,6)
        pade_ϵ3byf = get_pade(Polynomial(ϵ3byf),6,6)
        pade_ϵ4byf = get_pade(Polynomial(ϵ4byf),6,6)
        pade_byf = get_pade(Polynomial(byf),6,6)
        pade_gbyf2 = get_pade(Polynomial(gbyf),5,5)
        pade_ϵ2byf2 = get_pade(Polynomial(ϵ2byf),5,5)
        pade_ϵ3byf2 = get_pade(Polynomial(ϵ3byf),5,5)
        pade_ϵ4byf2 = get_pade(Polynomial(ϵ4byf),5,5)
        pade_byf2 = get_pade(Polynomial(byf),5,5)
        
        line_gbyf = CairoMakie.lines!(axchain,xvec, x->2*pade_gbyf(tanh(-f1*x)), color = (color_vec[1],0.7))
        line_ϵ2byf = CairoMakie.lines!(axchain,xvec, x->2*abs(pade_ϵ2byf(tanh(-f3*x))), color = (color_vec[2],0.7))
        line_ϵ3byf = CairoMakie.lines!(axchain,xvec, x->2*abs(pade_ϵ3byf(tanh(-f3*x))), color = (color_vec[4],0.7))
        line_ϵ4byf = CairoMakie.lines!(axchain,xvec, x->2*abs(pade_ϵ4byf(tanh(-f3*x))), color = (color_vec[5],0.7))
        line_byf = CairoMakie.lines!(axchain,xvec, x->2*abs(pade_byf(tanh(-f2*x))), color = (color_vec[3],0.7))
        CairoMakie.lines!(axchain,xvec, x->2*pade_gbyf2(tanh(-f1*x)), color = color_vec[1], linestyle = :dash)
        CairoMakie.lines!(axchain,xvec, x->2*abs(pade_ϵ2byf2(tanh(-f3*x))), color = color_vec[2], linestyle = :dash)
        CairoMakie.lines!(axchain,xvec, x->2*abs(pade_ϵ3byf2(tanh(-f3*x))), color = color_vec[4], linestyle = :dash)
        CairoMakie.lines!(axchain,xvec, x->2*abs(pade_ϵ4byf2(tanh(-f3*x))), color = color_vec[5], linestyle = :dash)
        CairoMakie.lines!(axchain,xvec, x->2*abs(pade_byf2(tanh(-f2*x))), color = color_vec[3], linestyle = :dash)

        CairoMakie.scatter!(axchain,QMCbetas, 2*gs./fs , color = color_vec[1])
        CairoMakie.scatter!(axchain,QMCbetas, 2*ϵ3s./fs, color = color_vec[4])
        CairoMakie.scatter!(axchain,QMCbetas, 2*ϵ4s./fs, color = color_vec[5])
        CairoMakie.scatter!(axchain,QMCbetas, 2*ϵ2s./fs, color = color_vec[2])
        CairoMakie.scatter!(axchain,QMCbetas,2* QMCbetas.*ones(length(fs))./fs, color = color_vec[3])
    end

    # Triangular lattice plot
    if true
        f1 = 0.4
        f2 = 0.4
        f3 = 0.4
        transf1 = get_LinearTrafoToCoeffs_u(12,f1)
        transf2 = get_LinearTrafoToCoeffs_u(13,f2)
        transf3 = get_LinearTrafoToCoeffs_u(12,f3)

        dists, indices = get_unique_distances(triang_lattice)
        invcorrs = invcorrs_triang
        gbyf = transf1*(invcorrs[indices[2][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ2byf = transf3*(invcorrs[indices[3][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ3byf = transf3*(invcorrs[indices[4][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ4byf = transf3*(invcorrs[indices[5][1]]/invcorrs[indices[1][1]]).coeffs
        byf = transf2*vcat(0,(1/invcorrs[indices[1][1]]).coeffs)[:]

        pade_gbyf = get_pade(Polynomial(gbyf),6,6)
        pade_ϵ2byf = get_pade(Polynomial(ϵ2byf),6,6)
        pade_ϵ3byf = get_pade(Polynomial(ϵ3byf),6,6)
        pade_ϵ4byf = get_pade(Polynomial(ϵ4byf),6,6)
        pade_byf = get_pade(Polynomial(byf),6,6)
        pade_gbyf2 = get_pade(Polynomial(gbyf),5,5)
        pade_ϵ2byf2 = get_pade(Polynomial(ϵ2byf),5,5)
        pade_ϵ3byf2 = get_pade(Polynomial(ϵ3byf),4,4)
        pade_ϵ4byf2 = get_pade(Polynomial(ϵ4byf),5,5)
        pade_byf2 = get_pade(Polynomial(byf),5,5)

        CairoMakie.lines!(axtriang,xvec, x->3*pade_gbyf(tanh(-f1*x)), color = (color_vec[1],0.7))
        CairoMakie.lines!(axtriang,xvec, x->3*abs(pade_ϵ2byf(tanh(-f3*x))), color = (color_vec[2],0.7))
        line_ϵ3byf = CairoMakie.lines!(axtriang,xvec, x->3*abs(pade_ϵ3byf(tanh(-f3*x))), color = (color_vec[4],0.7))
        line_ϵ4byf =CairoMakie.lines!(axtriang,xvec, x->3*abs(pade_ϵ4byf(tanh(-f3*x))), color = (color_vec[5],0.7))
        CairoMakie.lines!(axtriang,xvec, x->3*abs(pade_byf(tanh(-f2*x))), color = (color_vec[3],0.7))
        CairoMakie.lines!(axtriang,xvec, x->3*pade_gbyf2(tanh(-f1*x)), color = color_vec[1], linestyle = :dash)
        CairoMakie.lines!(axtriang,xvec, x->3*abs(pade_ϵ2byf2(tanh(-f3*x))), color = color_vec[2], linestyle = :dash)
        CairoMakie.lines!(axtriang,xvec, x->3*abs(pade_ϵ3byf2(tanh(-f3*x))), color = color_vec[4], linestyle = :dash)
        CairoMakie.lines!(axtriang,xvec, x->3*abs(pade_ϵ4byf2(tanh(-f3*x))), color = color_vec[5], linestyle = :dash)
        CairoMakie.lines!(axtriang,xvec, x->3*abs(pade_byf2(tanh(-f2*x))), color = color_vec[3], linestyle = :dash)

        legline = CairoMakie.lines!(axtriang,[-1,-2],[1,1], color = :black)
        legdashline = CairoMakie.lines!(axtriang,[-1,-2],[1,1], color = :black, linestyle = :dash)
        legscatter = CairoMakie.scatter!(axtriang,[-1,-2],[1,1], color = :black)
    end

    # Kagome lattice plot
    if true
        f1 = 0.8
        f2 = 0.8
        f3 = 0.8
        transf1 = get_LinearTrafoToCoeffs_u(12,f1)
        transf2 = get_LinearTrafoToCoeffs_u(13,f2)
        transf3 = get_LinearTrafoToCoeffs_u(12,f3)

        dists, indices = get_unique_distances(kagome_lattice)
        invcorrs = invcorrs_kagome
        gbyf = transf1*(invcorrs[indices[2][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ2byf = transf3*(invcorrs[indices[3][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ3byf = transf3*(invcorrs[indices[4][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ4byf = transf3*(invcorrs[indices[5][1]]/invcorrs[indices[1][1]]).coeffs
        byf = transf2*vcat(0,(1/invcorrs[indices[1][1]]).coeffs)[:]
        
        pade_gbyf = get_pade(Polynomial(gbyf),6,6)
        pade_ϵ2byf = get_pade(Polynomial(ϵ2byf),6,6)
        pade_ϵ3byf = get_pade(Polynomial(ϵ3byf),6,6)
        pade_ϵ4byf = get_pade(Polynomial(ϵ4byf),6,6)
        pade_byf = get_pade(Polynomial(byf),6,6)
        pade_gbyf2 = get_pade(Polynomial(gbyf),5,5)
        pade_ϵ2byf2 = get_pade(Polynomial(ϵ2byf),5,5)
        pade_ϵ3byf2 = get_pade(Polynomial(ϵ3byf),5,5)
        pade_ϵ4byf2 = get_pade(Polynomial(ϵ4byf),5,5)
        pade_byf2 = get_pade(Polynomial(byf),5,5)
        
        CairoMakie.lines!(axkagome,xvec, x->2*pade_gbyf(tanh(-f1*x)), color = (color_vec[1],0.7))
        CairoMakie.lines!(axkagome,xvec, x->2*abs(pade_ϵ2byf(tanh(-f3*x))), color = (color_vec[2],0.7))
        CairoMakie.lines!(axkagome,xvec, x->2*abs(pade_ϵ3byf(tanh(-f3*x))), color = (color_vec[4],0.7))
        CairoMakie.lines!(axkagome,xvec, x->2*abs(pade_ϵ4byf(tanh(-f3*x))), color = (color_vec[5],0.7))
        CairoMakie.lines!(axkagome,xvec, x->2*abs(pade_byf(tanh(-f2*x))), color = (color_vec[3],0.7))
        CairoMakie.lines!(axkagome,xvec, x->2*pade_gbyf2(tanh(-f1*x)), color = color_vec[1], linestyle = :dash)
        CairoMakie.lines!(axkagome,xvec, x->2*abs(pade_ϵ2byf2(tanh(-f3*x))), color = color_vec[2], linestyle = :dash)
        CairoMakie.lines!(axkagome,xvec, x->2*abs(pade_ϵ3byf2(tanh(-f3*x))), color = color_vec[4], linestyle = :dash)
        CairoMakie.lines!(axkagome,xvec, x->2*abs(pade_ϵ4byf2(tanh(-f3*x))), color = color_vec[5], linestyle = :dash)
        CairoMakie.lines!(axkagome,xvec, x->2*abs(pade_byf2(tanh(-f2*x))), color = color_vec[3], linestyle = :dash)
    end

    # Pyrochlore lattice plot
    if true
        f1 = 0.4
        f2 = 0.4
        f3 = 0.4
        transf1 = get_LinearTrafoToCoeffs_u(12,f1)
        transf2 = get_LinearTrafoToCoeffs_u(13,f2)
        transf3 = get_LinearTrafoToCoeffs_u(12,f3)
        dists, indices = get_unique_distances(pyrochlore_lattice)
        invcorrs = invcorrs_pyrochlore
        gbyf = transf1*(invcorrs[indices[2][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ2byf = transf3*(invcorrs[indices[3][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ3byf = transf3*(invcorrs[indices[4][1]]/invcorrs[indices[1][1]]).coeffs
        ϵ4byf = transf3*(invcorrs[indices[5][1]]/invcorrs[indices[1][1]]).coeffs
        byf = transf2*vcat(0,(1/invcorrs[indices[1][1]]).coeffs)[:]
        
        pade_gbyf = get_pade(Polynomial(gbyf),6,6)
        pade_ϵ2byf = get_pade(Polynomial(ϵ2byf),6,5)
        pade_ϵ3byf = get_pade(Polynomial(ϵ3byf),6,5)
        pade_ϵ4byf = get_pade(Polynomial(ϵ4byf),6,5)
        pade_byf = get_pade(Polynomial(byf),6,6)
        pade_gbyf2 = get_pade(Polynomial(gbyf),5,5)
        pade_ϵ2byf2 = get_pade(Polynomial(ϵ2byf),5,5)
        pade_ϵ3byf2 = get_pade(Polynomial(ϵ3byf),5,5)
        pade_ϵ4byf2 = get_pade(Polynomial(ϵ4byf),5,5)
        pade_byf2 = get_pade(Polynomial(byf),5,5)
        
        CairoMakie.lines!(axpyro,xvec, x->2*pade_gbyf(tanh(-f1*x)), color = (color_vec[1],0.7))
        CairoMakie.lines!(axpyro,xvec, x->2*abs(pade_ϵ2byf(tanh(-f3*x))), color = (color_vec[2],0.7))
        CairoMakie.lines!(axpyro,xvec, x->2*abs(pade_ϵ3byf(tanh(-f3*x))), color = (color_vec[4],0.7))
        CairoMakie.lines!(axpyro,xvec, x->2*abs(pade_ϵ4byf(tanh(-f3*x))), color = (color_vec[5],0.7))
        CairoMakie.lines!(axpyro,xvec, x->2*abs(pade_byf(tanh(-f2*x))), color = (color_vec[3],0.7))
        CairoMakie.lines!(axpyro,xvec, x->2*pade_gbyf2(tanh(-f1*x)), color = color_vec[1], linestyle = :dash)
        CairoMakie.lines!(axpyro,xvec, x->2*abs(pade_ϵ2byf2(tanh(-f3*x))), color = color_vec[2], linestyle = :dash)
        CairoMakie.lines!(axpyro,xvec, x->2*abs(pade_ϵ3byf2(tanh(-f3*x))), color = color_vec[4], linestyle = :dash)
        CairoMakie.lines!(axpyro,xvec, x->2*abs(pade_ϵ4byf2(tanh(-f3*x))), color = color_vec[5], linestyle = :dash)
        CairoMakie.lines!(axpyro,xvec, x->2*abs(pade_byf2(tanh(-f2*x))), color = color_vec[3], linestyle = :dash)
    end

    # Figure formatting
    hideydecorations!(axtriang)
    hideydecorations!(axpyro)
    hidexdecorations!(axchain)
    hidexdecorations!(axtriang)
    axtriang.ygridvisible =true
    axpyro.ygridvisible =true
    colgap!(grid, 4)
    rowgap!(grid, 2)

    # Legends
    leg = axislegend(axtriang, 
    [line_byf,line_gbyf, line_ϵ2byf, line_ϵ3byf, line_ϵ4byf ],
    [L"|x/f_x|",L"|g_x/f_x|", L"|\epsilon_2/f_x|", L"|\epsilon_3/f_x|", L"|\epsilon_4/f_x|"],
    tellheight = false,
    tellwidth = false,
    margin = (4, -1, 4, 4),
    framevisible = false,
    position = :rb,
    linewidth = 5,
    labelsize = 12,
    nbanks = 2, 
    rowgap = 2,
    colgap = 10,
    patchsize=(15,10))

    leg = axislegend(axkagome, 
    [legline, legdashline, legscatter],
    ["u-pade [6,6]","u-pade [5,5]","QMC"],
    tellheight = false,
    tellwidth = false,
    margin = (4, -1, 4, 4),
    framevisible = false,
    position = :rb,
    linewidth = 5,
    labelsize = 10,
    rowgap = -4)

    resize_to_layout!(fig)

    # Labels
    CairoMakie.text!(axchain, 0.2, 1.0, text=L"\times 2, \; f = 0.4", fontsize=11,color=:black, align=(:left, :top) )
    CairoMakie.text!(axtriang, 0.2, 1.0, text=L"\times 3, \; f = 0.4", fontsize=11,color=:black, align=(:left, :top) )
    CairoMakie.text!(axkagome, 0.2, 1.0, text=L"\times 2, \; f = 0.8", fontsize=11,color=:black, align=(:left, :top) )
    CairoMakie.text!(axpyro, 0.2, 1.0, text=L"\times 2, \; f = 0.4", fontsize=11,color=:black, align=(:left, :top) )

    # Save figure
    display(fig)
   # save("examples/Renormalized_Mean_Field_Theory/rMF.png",fig;px_per_unit = 600/96)
end