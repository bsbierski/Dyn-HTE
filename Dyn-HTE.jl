using JLD2
#using RobustPade


import Pkg 
Pkg.activate(@__DIR__) #activates the environment in the folder of the current file

include("Embedding.jl")
include("GraphGeneration.jl")
include("LatticeGraphs.jl")
include("ConvenienceFunctions.jl") 
#specify max order
max_order = 4

#LOAD FILES 
#-------------------------------------------------------------------------------------

#generate list of graphs
graphs_vec = [load_object("GraphFiles/graphs_"*string(nn)*".jld2") for nn in 0:max_order];
gG_vec = getGraphsG(graphs_vec);
## identify which gG have the same underlying simple-graph structure. Precalculate Symmetry factors. 
gG_vec_unique = give_unique_gG_vec(gG_vec);

#create vector of all lower order dictionaries
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,max_order+1) 
   
#load dictionaries of all lower orders C_Dict_vec 
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object("GraphEvaluations/Spin_S1half/C_"*string(ord)*".jld2")
end 
#-----------------------------------

#1. Define lattice ball for embedding (it is enough for embedding of max_order graphs to have ball radius L=max_order)
L = 4
lattice,LatGraph,center_sites = getLattice_Ball(L,"kagome");
display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.1,fontsize=7,nodeshape=:rect,curves=false))

#2.Compute all correlations in the lattice
Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);

#3. Fourier Transform



### Compute A 2D Brillouin zone cut: 
N = 100
kx = range(-8pi,8pi,length=N)
ky = range(-8pi,8pi,length=N) #[0.] #for chains
kmat = [(y,x,x) for x in kx, y in ky ]
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);

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
path = [(0,0),(2pi,0)]
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

x = 2.0
#k_vec = [(k,0.0) for k in 0.01:(2*π-0.02)/14:(2*π-0.01)]
w_vec = collect(0.01:0.0314/2:5.0)
JSkw_mat = get_JSkw_mat_finitex("pade",x,kvec,w_vec,0.02,0,3,200,false,Correlators,lattice,center_sites)





fig = Figure(fontsize=25,resolution = (900,500));
ax=Axis(fig[1,1],limits=(0,Nk+1,0,5),ylabel=L"\omega/J=w",title="Square Lattice: x="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25);
hm=CairoMakie.heatmap!(ax,[k for k in 1:Nk+1],w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,1.0),highclip=:white);
ax.xticks = (kticks_positioins,pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=40, label = L"J S(k,w)")
resize_to_layout!(fig);
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