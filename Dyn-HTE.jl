using JLD2
#using RobustPade


import Pkg 
Pkg.activate(@__DIR__) #activates the environment in the folder of the current file

include("Embedding.jl")
include("LatticeGraphs.jl")
include("ConvenienceFunctions.jl") 

#1. Define lattice ball for embedding (it is enough for embedding of max_order graphs to have ball radius L=max_order)
L = 10
spin_length = 1/2
hte_graphs = load_dyn_hte_graphs(spin_length,L);
hte_lattice = getLattice(L,"simple_cubic");
display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

@time c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs; verbose =false, max_order = 12);

#3. Fourier Transform

### Compute A 2D Brillouin zone cut: 
N = 50
kx = range(-2pi,2pi,length=N)
ky = range(-2pi,2pi,length=N) #[0.] #for chains
kmat = [(y,x,0.) for x in kx, y in ky ];
c_kDyn_mat =  get_c_kDyn(kmat,c_iipDyn_mat,hte_lattice);

x = -1.0
padetype = [5,5]
evaluate(y) = eval_correlator_LR_continuous_pad(y, x, padetype); #define evaluation function
struc = real.(evaluate.(c_kDyn_mat));
p = Plots.heatmap(kx,ky,struc)



### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
iomega = 0
JoverT = -10.5
padetype = [4,4]
evaluate(x) = eval_correlator_LR_continuous_pad_Mats(x,iomega, JoverT, padetype); #define evaluation function
struc = real.(evaluate.(structurefactor));
p = Plots.heatmap(kx,ky,struc,clims=(0,1.))

############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone

#---chain
path = [(0.01,0.01),(2pi-0.1,0.01)]
pathticks = ["Γ","Γ"]

#---triangular
path = [(0,0),(0,2pi/(sqrt(3))),(2pi/(3),2pi/(sqrt(3))),(0.,0.)]
pathticks = ["Γ","M","K","Γ"]

#---triangular Sherman 1
path = [(0,2pi/(sqrt(3))),(0.001,0.001),(4pi/(3),0),(1.0*pi,pi/sqrt(3)),(2pi/(3),0),(pi/2,pi/(2*sqrt(3)))]
pathticks = ["M","Γ","K","M","Y1","Y"]

#---triangular Sherman 2
path = [(0.1,0),(4pi/(3),0),(1.0*pi,pi/sqrt(3)),(4pi/(3),0),(0.1,0)]
pathticks = ["Γ","K","M","K","Γ"]

#---square
path = [(0,0),(pi,0),(pi,pi),(0,0)]
pathticks = ["Γ","X","M","Γ"]

#---square path Sherman
path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
pathticks = ["K","X","M","K","Γ","X"]


#Generate the path 
Nk = 100
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);

kticks_positioins
kvec[63]
BrillPath

####Plot the susceptibility along the high symmetry path for various temperatures 
#= Tvec = [2.0,1.0,0.5,0.375]
padetype = [4,4]

p= plot();
for T in Tvec
evaluate(x) = eval_correlator_LR_continuous_pad_Mats(x,0,-1/T,padetype);
cut = 1/T* real.(evaluate.(BrillPath));
plot!(p,1:length(kvec),cut, label = "T = $T", xticks=(kticks_positioins,pathticks)#= ,yrange = [0,2],xrange = [900,2800] =#)
end
p

BrillPath[1]

get_intDiffApprox(p::Polynomial,Tvec,4,2,2)


####Plot the Spectrum along the high symmetry path for a set temperature  
ωvec = -6:0.1:6
JoverT = 0.5
padetype = [4,4]
spectrum = transpose([imag(eval_correlator_LR_continuous_pad(x,ω,JoverT,padetype)) for x in BrillPath, ω in ωvec]);
heatmap(1:length(kvec),ωvec,spectrum, xticks=(kticks_positioins,pathticks),colormap = :RdBu )

 =#

###### S(k,w) heatmap}
using CairoMakie

x = 0
#k_vec = [(k,0.0) for k in 0.01:(2*π-0.02)/14:(2*π-0.01)]
w_vec = collect(0.01:0.0314/2:5.0)
JSkw_mat = get_JSkw_mat_finitex("total","pade",x,kvec,w_vec,0.02,1,5,200,false,Correlators,lattice,center_sites)






fig = Figure(fontsize=25,resolution = (900,500));
ax=Axis(fig[1,1],limits=(0,Nk+1,0,5),ylabel=L"\omega/J=w",title="Square Lattice: x="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25);
hm=CairoMakie.heatmap!(ax,[k for k in 1:Nk+1],w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,1.0),highclip=:white);
ax.xticks = (kticks_positioins,pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=40, label = L"J S(k,w)")
resize_to_layout!(fig);
display(fig)

CairoMakie.plot!(ax,QMC_axis*(Nk+0.5),1*QMCdata, color = :pink, alpha = 0.45,label = L"T=0 \text{QMC-Gap}")
CairoMakie.plot!(ax,[k for k in 1:Nk+1],2.4*disp, color = :orange, alpha = 0.45,label = L"T=0 \text{LSWT-Dispersion}")
axislegend(ax)
display(fig)

save("SquareX"*string(x)*".pdf",fig)


#plot(JSkw_mat[40,:])



display(fig)
kvec
disp = [sqrt(1-1/4*(cos(k[1])+cos(k[2]))^2)  for k in kvec]

I

using LinearAlgebra

# Function to compute matrix-valued Taylor coefficients for exp(xA)
function matrix_exponential_taylor(A, order)
    C =  Matrix{Float64}[]
    push!(C, I(2))  # C_0 = I
    for k in 1:order
        push!(C, A^k / factorial(k))  # C_k = A^k / k!
    end
    return C
end



using LinearAlgebra

# Function to compute the matrix-valued Padé approximant
function matrix_pade(taylor_coeffs, m, n, z)
    # taylor_coeffs is an array of matrices [A_0, A_1, A_2, ..., A_{m+n}]
    # m is the degree of the numerator polynomial
    # n is the degree of the denominator polynomial
    # z is the point at which to evaluate the Padé approximant

    # Construct the system of equations to solve for the coefficients of Q_n(z)
    # We use the fact that the Taylor series of R_{m,n}(z) should match the given Taylor series up to order m+n

    # Number of coefficients in the Taylor series
    L = length(taylor_coeffs)

    # Check if we have enough coefficients
    if L < m + n + 1
        error("Not enough Taylor coefficients to compute the Padé approximant.")
    end

    # Construct the matrix for the linear system
    # The system is of the form: C * q = b
    # where q is the vector of coefficients of Q_n(z), and b is derived from the Taylor coefficients

    # Initialize the matrix C and vector b
    C = zeros(eltype(taylor_coeffs[1]), (m+1)*size(taylor_coeffs[1], 1), (n+1)*size(taylor_coeffs[1], 2))
    b = zeros(eltype(taylor_coeffs[1]), (m+1)*size(taylor_coeffs[1], 1), size(taylor_coeffs[1], 2))

    # Fill the matrix C and vector b
    for i in 1:m+1
        for j in 1:n+1
            if i >= j
                C[(i-1)*size(taylor_coeffs[1], 1)+1:i*size(taylor_coeffs[1], 1), (j-1)*size(taylor_coeffs[1], 2)+1:j*size(taylor_coeffs[1], 2)] = taylor_coeffs[i-j+1]
            end
        end
        b[(i-1)*size(taylor_coeffs[1], 1)+1:i*size(taylor_coeffs[1], 1), :] = taylor_coeffs[i]
    end

    # Solve the linear system C * q = b
    q = C \ b

    # Extract the coefficients of Q_n(z)
    Q_coeffs = [reshape(q[(i-1)*size(taylor_coeffs[1], 2)+1:i*size(taylor_coeffs[1], 2), :], size(taylor_coeffs[1])) for i in 1:n+1]

    # Compute the coefficients of P_m(z)
    P_coeffs = [sum(taylor_coeffs[k+1] * Q_coeffs[j+1] for j in 0:min(k, n)) for k in 0:m]
    print(P_coeffs)
    print(Q_coeffs)
    # Evaluate P_m(z)
    P_z = sum(z^k * P_coeffs[k+1] for k in 0:m)

    # Evaluate Q_n(z)
    Q_z = sum(z^k * Q_coeffs[k+1] for k in 0:n)

    # Compute the Padé approximant R_{m,n}(z) = P_m(z) * inv(Q_n(z))
    R_z = P_z * inv(Q_z)

    return R_z
end

# Test the function with an easy example
# Let's consider a simple Taylor series with matrix coefficients: A_0 = I, A_1 = A, A_2 = A^2/2, etc.
A = [1 0;0 1]
taylor_coeffs = matrix_exponential_taylor(A,4)

# Choose m = 2, n = 2
m = 2
n = 2

# Evaluate at z = 1.0
z = 0.9

# Compute the Padé approximant


R_z = matrix_pade(taylor_coeffs, 2, 2, z)



println("Padé approximant at z = $z:")
println(R_z)


function scalar_pade(taylor_coeffs, m, n, z)
    # taylor_coeffs is an array of scalar coefficients [a_0, a_1, a_2, ..., a_{m+n}]
    # m is the degree of the numerator polynomial
    # n is the degree of the denominator polynomial
    # z is the point at which to evaluate the Padé approximant

    L = length(taylor_coeffs)
    if L < m + n + 1
        error("Not enough Taylor coefficients to compute the Padé approximant.")
    end

    # Construct the matrix for the linear system C * q = b
    C = zeros(m+1, n+1)
    b = zeros(m+1)

    for i in 1:m+1
        for j in 1:n+1
            if i >= j
                C[i, j] = taylor_coeffs[i-j+1]
            end
        end
        b[i] = taylor_coeffs[i]
    end

    # Solve the linear system C * q = b
    q = C \ b

    # Extract the coefficients of Q_n(z)
    Q_coeffs = q

    # Compute the coefficients of P_m(z)
    P_coeffs = [sum(taylor_coeffs[k+1] * Q_coeffs[j+1] for j in 0:min(k, n)) for k in 0:m]

    # Evaluate P_m(z)
    P_z = sum(z^k * P_coeffs[k+1] for k in 0:m)

    # Evaluate Q_n(z)
    Q_z = sum(z^k * Q_coeffs[k+1] for k in 0:n)

    # Compute the Padé approximant R_{m,n}(z) = P_m(z) / Q_n(z)
    R_z = P_z / Q_z

    return R_z
end



# Define a diagonal matrix A
A = Diagonal([1.0, 2.0])

# Taylor series coefficients for exp(A)
taylor_coeffs = [A^k / factorial(k) for k in 0:10]  # Up to order 4

# Choose m = 2, n = 2
m = 4
n = 4

# Evaluate at z = 1.0
z = 2.0

# Compute the matrix-valued Padé approximant
R_z_matrix = matrix_pade(taylor_coeffs, m, n, z)

println("Matrix-valued Padé approximant at z = $z:")
println(R_z_matrix)

# Compute the scalar Padé approximant for each diagonal element
R_z_scalar_1 = scalar_pade([1.0, 1.0, 1.0/2, 1.0/6, 1.0/24], m, n, z)  # For the first diagonal element (1.0)
R_z_scalar_2 = scalar_pade([1.0, 2.0, 4.0/2, 8.0/6, 16.0/24], m, n, z)  # For the second diagonal element (2.0)

println("Scalar Padé approximant for the first diagonal element (1.0):")
println(R_z_scalar_1)

println("Scalar Padé approximant for the second diagonal element (2.0):")
println(R_z_scalar_2)

# Compare the results
println("Comparison:")
println("Matrix-valued Padé approximant diagonal elements: ", diag(R_z_matrix))
println("Scalar Padé approximants: ", [R_z_scalar_1, R_z_scalar_2])





x_QMC_vec = collect(1:8)
imax=20
m_vec = [0,1,2]
TGiip_m = zeros((length(x_QMC_vec),imax+1,length(m_vec)))
TGiip_m_err = zeros((length(x_QMC_vec),imax+1,length(m_vec)))

using HDF5
for x_QMC_pos in eachindex(x_QMC_vec) 
    x_QMC=x_QMC_vec[x_QMC_pos]
    fid = h5open("/Users/Schneider.Benedikt/LRZ Sync+Share/_Home/UNI/PHD_Project/Dynamic High Temperature Series Expansion/QMC_Worm/SpinHalfAFMHeisenbergChain/job_BSb_beta"*string(x_QMC)*".out.h5", "r")
    for (m_pos,m) in enumerate(m_vec)
        TGiip_m[x_QMC_pos,:,m_pos] = read(fid["simulation"]["results"]["DensDens_CorrFun_w$m"]["mean"]["value"])[1:imax+1]
        TGiip_m_err[x_QMC_pos,:,m_pos] = read(fid["simulation"]["results"]["DensDens_CorrFun_w$m"]["mean"]["error"])[1:imax+1]
    end
end

gG_vec_unique_new = unique_Graphs(gG_vec_unique)
@save "gG_vec_unique_new.jld2" gG_vec_unique_new

#= 
function compute_lattice_correlations_new(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec)::Array{Matrix{Rational{Int64}}}
    """compute all correlations from the center_sites to all other sites of the lattice"""
    Correlators = Array{Matrix{Rational{Int64}}}(undef, lattice.length,length(lattice.unitcell.basis));
    Threads.@threads for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
        Correlators[jp,b] = mapreduce(permutedims, vcat, Calculate_Correlator_fast_new(LatGraph,center_sites[b],jp,max_order,gG_vec_unique,C_Dict_vec))
        end
    end
    return Correlators
end

@time Correlators_new = compute_lattice_correlations_new(LatGraph,lattice,center_sites,max_order,gG_vec_unique_new,C_Dict_vec)


for maxorder = 0:0
file_path = "GraphFiles/unique_gG_vec_$maxorder.jld2"
if isfile(file_path)
        unique_gG_vec = load_object(file_path) 
        unique_gG_vec = unique_Graphs(unique_gG_vec)
       @save "GraphFiles/unique_gG_vec_$maxorder"*".jld2" unique_gG_vec
end
end

print("test$test"*"test")

maxorder = 12
file_path = "GraphFiles/unique_gG_vec_$maxorder.jld2"
if isfile(file_path)
    unique_gG_vec = load_object(file_path) 
    gG_vec_unique_new = unique_Graphs(unique_gG_vec)
   #@save "GraphFiles/unique_gG_vec_$maxorder"*"_new.jld2"
end

@load file_path gG_vec_unique_new

gG_vec_unique_new.graphs

unique_Graphs(12,)
vcat([1],[2],[3],[4])



unique_graphs_12_part = unique_Graphs(12,gG_vec_unique.graphs[1:5000])
@save "GraphFiles/unique_gG_vec_$maxorder"*"_1"*".jld2" unique_graphs_12_part
unique_graphs_12_part = unique_Graphs(12,gG_vec_unique.graphs[5001:40000])
@save "GraphFiles/unique_gG_vec_$maxorder"*"_2"*".jld2" unique_graphs_12_part
unique_graphs_12_part = unique_Graphs(12,gG_vec_unique.graphs[40001:150000])
@save "GraphFiles/unique_gG_vec_$maxorder"*"_3"*".jld2" unique_graphs_12_part
unique_graphs_12_part = unique_Graphs(12,gG_vec_unique.graphs[150001:end]
@save "GraphFiles/unique_gG_vec_$maxorder"*"_4"*".jld2" unique_graphs_12_part
for part = 1:6
split = split_vec(gG_vec_unique.graphs,part,10)[1]
unique_graphs_12_part = unique_Graphs(12,split)
@save "GraphFiles/unique_gG_vec_$maxorder"*"_$part"*"_new.jld2" unique_graphs_12_part
end

split = split_vec(gG_vec_unique_new.graphs,part,4)[1]
unique_Graphs([12,split]) =#
max_order = 12
gG_vec_unique = give_unique_gG_vec(max_order);
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,max_order+1) ;
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object("GraphEvaluations/Spin_S1half/C_"*string(ord)*".jld2")
end 

@time gG_vec_unique_precalc = precalculate_unique_graphs(max_order,gG_vec_unique,C_Dict_vec);

gG_vec_unique_precalc.graphs[12].graph_value


function precalculate_unique_graphs(max_order::Int,gG_vec_unique::unique_Graphs,C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}})::unique_Graphs_precalc
    """ Calculate the coefficients of (-x)^n for TG_ii'(iν_m) from embedding factors of only the unique simple graphs and the gG's symmetry factors """    

    result_vector = Vector{unique_Graph_precalc}(undef,length(gG_vec_unique.graphs))

    # only iterate over the unique simple graphs in unique_Gg
    for (index,unique_Gg) in enumerate(gG_vec_unique.graphs)
            #initialize result array
        result_array = Matrix{Rational{Int128}}(undef, max_order+1,10)

     #for every order we get result vector representing prefactors of [δw,Δ^2,Δ^4,Δ^6,Δ^8,Δ^10,Δ^12,Δ^14,Δ^16,Δ^18]
     for ord = 1:max_order+1
        result_array[ord,:] = zeros(Rational{Int128},10)
    end

        gg = unique_Gg.ref_graph   #Graph
        gg_dist = unique_Gg.distance

        

        #### now we sum overall graphG eqivalent to the unique Gg
        for graph in unique_Gg.gG_vec
            g_order = graph.order #order
            gG_vec_index = graph.index #index
            symmetry_factor = graph.symmetry_factor#symmetry factor
            is_symmetric = graph.is_symmetric  #bool if the graph is symmetric
           
            fac = 2
            if is_symmetric
                fac = 1
            end

            #look up the value of the graph from C_Dict_vec
            look_up_dict =C_Dict_vec[g_order+1][gG_vec_index]
            result_array[g_order+1,:] .+= look_up_dict/symmetry_factor*fac
        end

        result_vector[index] = unique_Graph_precalc(gg,gg_dist,result_array)
    end

    return unique_Graphs_precalc(max_order,result_vector)
end


function compute_lattice_correlations_new(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec)::Array{Matrix{Rational{Int128}}}
    """compute all correlations from the center_sites to all other sites of the lattice"""
    Correlators = Array{Matrix{Rational{Int128}}}(undef, lattice.length,length(lattice.unitcell.basis));
    Threads.@threads for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
        Correlators[jp,b] =  Calculate_Correlator_faster(LatGraph,center_sites[b],jp,max_order,gG_vec_unique,C_Dict_vec)
        end
    end
    return Correlators
end

function Calculate_Correlator_faster(L::SimpleGraph{Int},ext_j1::Int,ext_j2::Int,max_order::Int,gG_vec_unique::unique_Graphs_precalc,C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}})::Matrix{Rational{Int128}}
    """ Calculate the coefficients of (-x)^n for TG_ii'(iν_m) from embedding factors of only the unique simple graphs and the gG's symmetry factors """    

    #initialize result array
    result_array = Matrix{Rational{Int128}}(undef, max_order+1,10)

    #for every order we get result vector representing prefactors of [δw,Δ^2,Δ^4,Δ^6,Δ^8,Δ^10,Δ^12,Δ^14,Δ^16,Δ^18]
    for ord = 1:max_order+1
        result_array[ord,:] = zeros(Rational{Int128},10)
    end

    #calculate the shortest graph distance between ext_j1 and ext_j2
    ext_dist = dijkstra_shortest_paths(L,ext_j1).dists[ext_j2]

    # only iterate over the unique simple graphs in unique_Gg
    for  unique_Gg in gG_vec_unique.graphs
        gg = unique_Gg.ref_graph   #Graph
        gg_dist = unique_Gg.distance #edge distance between the external vertices
        gg_val = unique_Gg.graph_value
          # if the graph is long enough
          if gg_dist < ext_dist 
            continue
        end

        #if its the onsite correlator we only need on-site graphs 
        if ext_dist == 0
            if gg_dist > ext_dist 
                continue
            end
        else #if not, we dont need any on site graphs
            if gg_dist == 0 
                continue
            end
        end

        #calculate the embedding factor
        emb_fac = e_fast(L,ext_j1,ext_j2,gg)

        result_array .+=  emb_fac*gg_val
    end

    return result_array
end


max_order = 10

#LOAD FILES 
#-------------------------------------------------------------------------------------
#load list of unique graphs
gG_vec_unique = give_unique_gG_vec(max_order);

#create vector of all lower order dictionaries
C_Dict_vec = Vector{Vector{Vector{Rational{Int128}}}}(undef,max_order+1) ;
#load dictionaries of all lower orders C_Dict_vec 
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object("GraphEvaluations/Spin_S1half/C_"*string(ord)*".jld2")
end 
#-----------------------------------

#1. Define lattice ball for embedding (it is enough for embedding of max_order graphs to have ball radius L=max_order)
L = max_order
lattice,LatGraph,center_sites = getLattice_Ball(L,"triang");
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.1,fontsize=7,nodeshape=:rect,curves=false))

#2.Compute all correlations in the lattice
@time Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);

@time gG_vec_unique_precalc = precalculate_unique_graphs(max_order,gG_vec_unique,C_Dict_vec);
@time Correlators_new = compute_lattice_correlations_new(LatGraph,lattice,center_sites,max_order,gG_vec_unique_precalc,C_Dict_vec);

Correlators == Correlators_new


gG_vec_unique_precalc.graphs[2].graph_value




#Calculate the Quantum to Classical Correspondence


function calc_taylorinvmat_fun(corr)
    orderJ,orderω = size(corr)
    @variables x
    taylormat = y -> sum(y[i,1]*x^(i-1) for i =  1:orderJ)
    invtaylormat = inv(taylormat(corr));
    t = Taylor1(Float64,orderJ-1)
    return  substitute(invtaylormat, Dict(x=>t)) 
end


### Compute A 2D Brillouin zone cut: 
N = 30;
kx = (1:1)*2pi/(1);
ky = (1:N)*2pi/(N); #[0.] #for chains
kmat = [(y,x) for x in kx, y in ky ];
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);
invstruc = calc_taylorinvmat_fun.(structurefactor);
invcorrs = brillouin_zone_cut_inverse(kmat,invstruc,lattice,center_sites);


test = brillouin_zone_cut_inverse(kmat,structurefactor,lattice,center_sites);

test[12]
Float64.(Correlators[12])

sumcoeffs = [sum(abs.(invcorrs[i][:])) for i in 1:length(invcorrs)]
distances = [norm(lattice.sitePositions[i] .- lattice.sitePositions[center_sites[1]]) for i =1:length(lattice)]
pchain = scatter(psquare,distances,sumcoeffs, yaxis = :log, title = L"$G_{ij}^{-1}(i\omega = 0) = \sum_n c_{ij,n} x^n$" , xlabel = "|i-j|" , ylabel = L"\sum_n |c_{ij,n}|", label = "Chain")

(invcorrs[11]/invcorrs[13])[:]*10^6

invcorrs[13][1]

distances

center_sites

