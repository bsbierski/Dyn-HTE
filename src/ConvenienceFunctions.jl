"""
    create_spin_string(s::Number) -> String

Convert a numerical spin length to a standardized string representation.

# Arguments
- `s`: Spin length that can be converted to a rational number

# Returns
- A string in the format "S{n}" for integer spins or "S{n}half" for half-integer spins
  where {n} is the numerator of the rationalized spin length

# Examples
```julia
create_spin_string(1)    # returns "S1"
create_spin_string(1.5)  # returns "S3half"
create_spin_string(1/2)  # returns "S1half"
```

# Throws
- `ErrorException` if the spin length cannot be represented as an integer or half-integer value

# Notes 
- Input is automatically converted to a rational number using `rationalize(1.0*s)`
- Only supports integer and half-integer spin values (denominator must be 1 or 2)
"""
function create_spin_string(s)
    spin_length = rationalize(1.0*s)
    denom = denominator(spin_length)
    is_half = ""
    if denom ==2
        is_half = "half"
    elseif denom != 1
        error("invalid spin_length: "*string(spin_length))
    end
    return "S"*string(numerator(spin_length))*is_half
end

"""
    get_c_iipDyn_mat(Graph, basis_positions::Vector{<:Int}, hte_graphs::Dyn_HTE_Graphs; 
                    verbose=false, max_order=12) -> Array{Matrix{Rational{Int128}}}
    get_c_iipDyn_mat(hte_lattice::Dyn_HTE_Lattice, hte_graphs::Dyn_HTE_Graphs; 
                    verbose=false, max_order=12) -> Array{Matrix{Rational{Int128}}}

Calculate the high-temperature series expansion of dynamic correlators G_ii' between basis sites
and all sites in the graph.

# Arguments
- `Graph`: The graph representing the lattice structure (first method)
- `basis_positions::Vector{<:Int}`: Vector of indices for basis sites (first method)
- `hte_lattice::Dyn_HTE_Lattice`: Object containing the lattice structure, graph, and basis positions (second method). See [`Dyn_HTE_Lattice`](@ref).
- `hte_graphs::Dyn_HTE_Graphs`: Container of precomputed HTE graphs and coefficients. See [`Dyn_HTE_Graphs`](@ref).
- `verbose::Bool=false`: When true, prints progress information
- `max_order::Int=12`: Maximum expansion order (will be capped by available data)

# Returns
- `Array{Matrix{Rational{Int128}}}`: An array where element `[jp,b]` contains the expansion
  coefficients for correlations between site `jp` and basis site `b`. 
The expansion coefficients themselves are written in a Matrix whose rows indicate orders of ``-x = -J/T``, and columns correspond to Matsubara frequency ``m=0`` (first column) or even powers of ``\\Delta = 1/(2\\pi m)`` (other columns).

# Notes
- The second method provides significant performance improvements by leveraging lattice symmetries
  to reduce the number of correlation calculations needed

# Examples
```julia
# Standard method: Calculate correlators for a square lattice to 12th order for spin 1/2
L = 12
n_max = 1*L
spin_length = 1/2
hte_graphs = load_dyn_hte_graphs(spin_length,L)
hte_lattice = getLattice(L,"square")
c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice.graph, hte_lattice.basis_positions, hte_graphs)

# Optimized method: Calculate correlators with symmetry optimization
L = 12
spin_length = 1/2
hte_graphs = load_dyn_hte_graphs(spin_length, L)
hte_lattice = getLattice(L, "square")
c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice, hte_graphs, verbose=true)
```
"""
function get_c_iipDyn_mat(Graph,basis_positions::Vector{<:Int},hte_graphs::Dyn_HTE_Graphs; verbose =false, max_order = 12)::Array{Matrix{Rational{Int128}}}
    ##preallocate output matrix 
    GiipDyn_mat = Array{Matrix{Rational{Int128}}}(undef, nv(Graph),length(basis_positions));

    unique_graphs = hte_graphs.unique_graphs
    c_dict = hte_graphs.c_dict
    
    # restricts maximum_order 
    max_order = min(max_order,unique_graphs.max_order)

    ##compute correlations
    Threads.@threads :dynamic for jp = 1:nv(Graph)
        if verbose == true
            println("bond "*string(jp)*" of "*string(nv(Graph)))
        end

        for b = eachindex(basis_positions)
            GiipDyn_mat[jp,b] = mapreduce(permutedims, vcat, Calculate_Correlator_fast(Graph,basis_positions[b],jp,max_order,unique_graphs,c_dict))
        end

    end

    return GiipDyn_mat
end

function get_c_iipDyn_mat(hte_lattice::Dyn_HTE_Lattice,hte_graphs::Dyn_HTE_Graphs; verbose =false, max_order = 12)::Array{Matrix{Rational{Int128}}}
    ##try to use lattice symmetries 
    try 
        sym_G,transl_G = getSymmetryGroup(hte_lattice.name) 
    catch err
        println("Symmetry for lattice "*hte_lattice.name*"not implemented, continue without using symmetries:" )
        return get_c_iipDyn_mat(hte_lattice.graph,hte_lattice.basis_positions,hte_graphs; verbose = verbose, max_order = max_order)
    end

    println("Calculating symmetry relations")
    sym_G,transl_G = getSymmetryGroup(hte_lattice.name) 

    lattice = hte_lattice.lattice
    #use symmetries to reduce necessary bond calculations. 
    reduction_dict,bond_vec_red,position_dict = sym_reduced_lattice(lattice,hte_lattice.basis_positions,sym_G,transl_G)
    
    println("Symmetry relations calculated")
    println("Calculating c_iipDyn_mat")

    ##preallocate output matrix 
    GiipDyn_mat = Array{Matrix{Rational{Int128}}}(undef, lattice.length,length(lattice.unitcell.basis));
    reduced_Giip = Vector{Matrix{Rational{Int128}}}(undef, length(bond_vec_red));

    unique_graphs = hte_graphs.unique_graphs
    c_dict = hte_graphs.c_dict
      
    # restricts maximum_order 
    max_order = min(max_order,unique_graphs.max_order)

  

    #calculate correlators
    Threads.@threads for i in reverse(eachindex(bond_vec_red))
        if verbose == true
            println("bond "*string(i)*" of "*string(length(bond_vec_red)))
        end

        reduced_Giip[i] = mapreduce(permutedims, vcat, Calculate_Correlator_fast(hte_lattice.graph,position_dict[i][1],position_dict[i][2],max_order,unique_graphs,c_dict))
    end
    
    ## fill the full correlation matrix according to the reduction_dict
    for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
        GiipDyn_mat[jp,b] = reduced_Giip[reduction_dict[CartesianIndex(b,jp)]]
        end
    end

    return GiipDyn_mat
end

"""
    get_c_iipEqualTime_mat(c_iipDyn_mat::Matrix{Matrix{Rational{Int128}}}) -> Array{Vector{Rational{Int128}}}

Perform frequency summation over real-space dynamic correlators to obtain equal-time correlators.

This function converts dynamic frequency-dependent correlation functions into equal-time
correlations by performing the appropriate frequency summation with predefined coefficients. This allows to compare with conventional HTE.

# Arguments
- `c_iipDyn_mat::Matrix{Matrix{Rational{Int128}}}`: Matrix of dynamic correlators where each element 
  contains a matrix of expansion coefficients. See [`get_c_iipDyn_mat`](@ref).

# Returns
- `Array{Vector{Rational{Int128}}}`: Matrix with the same dimensions as input, but where each element 
  is a vector of equal-time correlation coefficients

# Examples
```julia
# Convert dynamic correlators to equal-time correlators
c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice, hte_graphs)
c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat)

```

# Notes
- The frequency sum is performed using fixed rational coefficients for each frequency component
- The output has dimensions `[length(c_iipDyn_mat[:,1]), length(c_iipDyn_mat[1,:])]`
- Each entry in the returned matrix is a vector of length `max_order_plus1`
"""
function get_c_iipEqualTime_mat(c_iipDyn_mat::Matrix{Matrix{Rational{Int128}}})::Array{Vector{Rational{Int128}}}
   
    max_order_plus1 = size(c_iipDyn_mat[1,1])[1]
    c_iipEqualTime_mat = Array{Vector{Rational{Int128}}}(undef, length(c_iipDyn_mat[:,1]), length(c_iipDyn_mat[1,:]))
    for j in eachindex(c_iipDyn_mat[:,1])
        for b in eachindex(c_iipDyn_mat[1,:])
            c_iipEqualTime_mat[j,b] = [sum(c_iipDyn_mat[j,b][n,:] .* [1//1,1//12,1//720,1//30240,1//1209600,1//47900160,691//1307674368000,1//74724249600,3617//10670622842880000,43867//5109094217170944000]) for n in 1:max_order_plus1]
        end
    end
    return c_iipEqualTime_mat
end



###### bare series polynomial in Gii'(x,m) at Matsubara integer m truncated at n 

"""
    get_TGiip_Matsubara_xpoly(c_iipDyn_mat::Matrix{Matrix{Rational{Int128}}}, i::Int, ip::Int, m::Int) -> Polynomial{Float64}

Convert a dynamic correlator to a Matsubara frequency polynomial representation.

Extracts the time-ordered Matsubara correlator TGii'(iνm) as a polynomial in x = J/T for the 
spatial indices (i,ip) at Matsubara frequency index m.

# Arguments
- `c_iipDyn_mat::Matrix{Matrix{Rational{Int128}}}`: Matrix of dynamic correlator coefficients obtained from [`get_c_iipDyn_mat`](@ref)
- `i::Int`: First spatial index for the correlator
- `ip::Int`: Second spatial index for the correlator
- `m::Int`: Matsubara frequency index (m=0 is handled as a special case)

# Returns
- `Polynomial{Float64}`: Polynomial representation of the Matsubara correlator in variable x

# See Also
- [`flipEvenIndexEntries`](@ref): Function used internally to convert between different series expansions
"""
function get_TGiip_Matsubara_xpoly(c_iipDyn_mat::Matrix{Matrix{Rational{Int128}}},i::Int,ip::Int,m::Int)
    n_max = size(c_iipDyn_mat[1,1])[1]-1
    if m==0
        p_x = 1.0*Polynomial(flipEvenIndexEntries(c_iipDyn_mat[i,ip][:,1]))
    else
        coeffs_m = [sum([c_iipDyn_mat[i,ip][n+1,lhalf+1] * 1/(2*π*m)^(2*lhalf) for lhalf in 1:9]) for n in 0:n_max]
        p_x = 1.0*Polynomial(flipEvenIndexEntries(coeffs_m))
    end

    return p_x
end

"""
    flipEvenIndexEntries(v) -> Vector

Apply alternating signs to vector elements based on their position.

Transforms a vector [a,b,c,d,...] into [+a,-b,+c,-d,...] by multiplying elements
at even indices by -1. This changes an expansion in ''-x'' to an expansion in ''x''.

# Arguments
- `v`: Input vector of any numeric type

# Returns
- A vector of the same type and length as the input, but with alternating signs

# Usage
This function is used throughout the module to convert between different series expansions:

# See Also
- [`get_TGiip_Matsubara_xpoly`](@ref): Uses this to process correlation matrix coefficients
- [`get_moments_from_c_kDyn`](@ref): Applies sign flipping when extracting moments from dynamic coefficients
"""
function flipEvenIndexEntries(v)
    signs = [-1*(-1)^n for n in eachindex(v)]
    return v .* signs
end

###### resummation tools for polynomial p

"""
    get_pade(p::Polynomial, N::Int, M::Int) -> Polynomial

Create a Padé approximant for polynomial p with specified degrees.

Constructs a robust [N,M] Padé approximant for the input polynomial, handling potential
numerical issues like division by zero.

# Arguments
- `p::Polynomial`: Input polynomial to approximate
- `N::Int`: Degree of the numerator polynomial
- `M::Int`: Degree of the denominator polynomial

# Returns
- Padé approximant as a rational function approximating the original function

# Examples
```julia
# Create and evaluate Padé approximant
x_vec = collect(0:0.05:5.0)
coeffs_x = flipEvenIndexEntries(get_c_k(k, c_iipEqualTime_mat, hte_lattice))
p_x = Polynomial(coeffs_x)

# Create different Padé approximants
pade_66 = get_pade(p_x, 6, 6)
pade_55 = get_pade(p_x, 5, 5)

# Evaluate at multiple points
y_vec_66 = pade_66.(x_vec)
y_vec_55 = pade_55.(x_vec)
```

# See Also
- [`flipEvenIndexEntries`](@ref): Used to convert between different series expansions
- [`get_c_k`](@ref): Fourier transform function used in the example
- [`get_c_iipEqualTime_mat`](@ref): Function to obtain equal-time correlators
"""
function get_pade(p::Polynomial,N::Int,M::Int)
    
    return robustpade(p,N,M)
end


###### variable transform from x to u=tanh(fx)
"""
    get_p_u(coeffs_x::Vector{Float64}, f::Float64) -> Polynomial

Transform a polynomial in variable x to a polynomial in u=tanh(fx).

This function performs a variable substitution from x to u=tanh(fx) in a polynomial
defined by its coefficients, truncating the result to the same degree as the input.

# Arguments
- `coeffs_x::Vector{Float64}`: Coefficients of the polynomial in variable x
- `f::Float64`: Scaling factor in the tanh transformation

# Returns
- `Polynomial`: Transformed polynomial in variable u=tanh(fx)

# Examples
```julia
# Transform a polynomial in x to a polynomial in u=tanh(fx)
x_poly = Polynomial([1.0, 2.0, 3.0])  # 1 + 2x + 3x²
f = 0.48
u_poly = get_p_u(coeffs(x_poly), f)

# For larger polynomials, use the more efficient linear transformation
# This is equivalent but much faster for larger degree polynomials:
ufromx_mat = get_LinearTrafoToCoeffs_u(length(coeffs(x_poly))-1, f)
u_poly_fast = Polynomial(ufromx_mat * coeffs(x_poly))
```

# Notes
- This implementation is relatively slow for high-degree polynomials
- For better performance with multiple transformations, use [`get_LinearTrafoToCoeffs_u`](@ref) instead
"""
function get_p_u(coeffs_x::Vector{Float64},f::Float64)
    @variables x u
    x = taylor(atanh(u)/f, u, 0:(length(coeffs_x)-1), rationalize=false)
    p_u_ext = simplify(series(coeffs_x,x);expand=true)
    p_u = Polynomial(Symbolics.value.(taylor_coeff(p_u_ext,u,0:12,rationalize=false)),:u)
    return p_u
end

"""
    get_LinearTrafoToCoeffs_u(max_order::Int, f::Float64) -> Matrix{Float64}

Get the linear transformation matrix to convert polynomial coefficients from x to u=tanh(fx).

Returns a matrix that, when multiplied with coefficients of a polynomial in x, 
gives the coefficients for the same polynomial expressed in u=tanh(fx).

# Arguments
- `max_order::Int`: Maximum order of the polynomial (must be ≤ 16)
- `f::Float64`: Scaling factor in the tanh transformation

# Returns
- `Matrix{Float64}`: Transformation matrix of size `(max_order+1) × (max_order+1)`

# Examples
```julia
# Transform a polynomial from x-space to u-space
f = 0.48
x_vec = 0:0.1:4.0
u_vec = tanh.(f .* x_vec)

# Get original polynomial coefficients
coeffs_x = [1.0, 0.5, 0.25]  # 1 + 0.5x + 0.25x²
p_x = Polynomial(coeffs_x)

# Transform to u-space
ufromx_mat = get_LinearTrafoToCoeffs_u(length(coeffs_x)-1, f)
p_u = Polynomial(ufromx_mat * coeffs_x)

# Compare evaluations
y_x = p_x.(x_vec)
y_u = p_u.(u_vec)  # Should match y_x
```

# Notes
- More efficient than [`get_p_u`](@ref) for repeated transformations
- Precomputed for polynomials up to order 16
- Used in u-Padé resummation techniques for better convergence
- The transformation is defined as x = atanh(u)/f

# Throws
- `ErrorException` if `max_order > 16`
"""
function get_LinearTrafoToCoeffs_u(max_order::Int, f::Float64)::Matrix{Float64}

    if max_order > 16
        throw(error("max_order must be smaller than 17"))
    end

    data = [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1/f, 0, 1/(3f), 0, 1/(5f), 0, 1/(7f), 0, 1/(9f), 0, 1/(11f), 0, 1/(13f), 0, 1/(15f), 0],
        [0, 0, 1/f^2, 0, 2/(3f^2), 0, 23/(45f^2), 0, 44/(105f^2), 0, 563/(1575f^2), 0, 3254/(10395f^2), 0, 88069/(315315f^2), 0, 11384/(45045f^2)],
        [0, 0, 0, 1/f^3, 0, 1/f^3, 0, 14/(15f^3), 0, 818/(945f^3), 0, 141/(175f^3), 0, 13063/(17325f^3), 0, 16774564/(23648625f^3), 0],
        [0, 0, 0, 0, 1/f^4, 0, 4/(3f^4), 0, 22/(15f^4), 0, 1436/(945f^4), 0, 21757/(14175f^4), 0, 11368/(7425f^4), 0, 35874836/(23648625f^4)],
        [0, 0, 0, 0, 0, 1/f^5, 0, 5/(3f^5), 0, 19/(9f^5), 0, 457/(189f^5), 0, 7474/(2835f^5), 0, 261502/(93555f^5), 0],
        [0, 0, 0, 0, 0, 0, 1/f^6, 0, 2/f^6, 0, 43/(15f^6), 0, 680/(189f^6), 0, 3982/(945f^6), 0, 147668/(31185f^6)],
        [0, 0, 0, 0, 0, 0, 0, 1/f^7, 0, 7/(3f^7), 0, 56/(15f^7), 0, 688/(135f^7), 0, 12926/(2025f^7), 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1/f^8, 0, 8/(3f^8), 0, 212/(45f^8), 0, 6568/(945f^8), 0, 18778/(2025f^8)],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/f^9, 0, 3/f^9, 0, 29/(5f^9), 0, 2897/(315f^9), 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/f^10, 0, 10/(3f^10), 0, 7/f^10, 0, 748/(63f^10)],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/f^11, 0, 11/(3f^11), 0, 374/(45f^11), 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/f^12, 0, 4/f^12, 0, 146/(15f^12)],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/f^13, 0, 13/(3f^13), 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/f^14, 0, 14/(3f^14)],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/f^15, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/f^16]
    ]

    mat = zeros(Float64, 17, 17)
    for i in 1:17
        for j in 1:length(data[i])
            mat[j, i] = data[i][j]
        end
    end

    return mat[1:max_order+1,1:max_order+1]
end

"""
    extrapolate_series(series, method::String, parameters) -> Polynomial

Extrapolate a polynomial using specified approximation method with given parameters.

This function provides a unified interface to different series extrapolation techniques,
applying the chosen method with provided parameters to the input series.

# Arguments
- `series`: Polynomial to extrapolate 
- `method::String`: Extrapolation method, either "pade" or "u_pade"
- `parameters`: Method-specific parameters:
  - For "pade": (N, M) degrees for numerator and denominator
  - For "u_pade": (N, M, f) where f is the scaling factor for tanh transformation

# Returns
- Extrapolated result as a rational function, either in x or in u

# Examples
```julia
# Extrapolate a Polynomial via a a [N,M] Padé approximant
poly = Polynomial([1,2,3,4,5,6])
extrapolated_poly = extrapolate_series(poly,"pade",[3,2]))

```

# Notes
- For "pade", applies a [N,M] Padé approximant directly to the series
- For "u_pade", first applies the substitution u = tanh(f·x) then constructs a Padé approximant
- The "u_pade" method usually is more stable for small temperatures, but trivially gets constant for even smaller temperatures due to the saturation of tanh(f·x)

# See Also
- [`get_pade`](@ref): Used for Padé approximant construction
- [`get_LinearTrafoToCoeffs_u`](@ref): Used for efficient variable transformation in "u_pade" method

"""
function extrapolate_series(series,method::String,parameters) 
    if method == "pade"
        return get_pade(series,Int64(parameters[1]),Int64(parameters[2]))
    elseif method == "u_pade"
        substitution_matrix = get_LinearTrafoToCoeffs_u(length(coeffs(series))-1,parameters[3])
        p_u = Polynomial(substitution_matrix*coeffs(series))
        return get_pade(p_u,Int64(parameters[1]),Int64(parameters[2]))
    end

end



###### k-space functions
"""
    create_brillouin_zone_path(points, num_samples::Int) -> (Vector{NTuple{N,Float64}}, Vector{Int})

Create a linear interpolation path between high-symmetry points in a Brillouin zone.

This function takes a list of points (typically high-symmetry points) and creates a path
through them with a specified number of sampling points, distributed proportionally to
segment lengths.

# Arguments
- `points`: List of points in the Brillouin zone, each as a tuple of coordinates
- `num_samples::Int`: Total number of points to generate along the path

# Returns
- `Vector{NTuple{N,Float64}}`: List of interpolated points along the path
- `Vector{Int}`: Indices of the original input points in the interpolated path

# Examples
```julia
# Create a path through high-symmetry points in a triangular lattice
Γ, K, M = (0,0), (2π/3, 2π/sqrt(3)), (0, 2π/sqrt(3))
path = [Γ, K, M, Γ]  # Path through high-symmetry points
pathticks = ["Γ", "K", "M", "Γ"]  # Labels for plotting

# Generate 200 points along this path
k_vec, kticks_positions = create_brillouin_zone_path(path, 200)

# Use for plotting band structure or dispersion relations
plt = plot(xlabel="k", ylabel="ω")
plot!(plt, 0:length(k_vec)-1, dispersion_function.(k_vec), 
      xticks=(kticks_positions .- 1, pathticks))
```

# Notes
- The number of points per segment is proportional to the segment length
- The original input points are always included in the output
- Useful for creating paths for spectral function plots or band structure plots

# See Also
- [`get_c_k`](@ref): Used to compute correlators along the generated k-path
- [`get_JSkw_mat`](@ref): Uses k-point paths for dynamic structure factor calculations
"""
function create_brillouin_zone_path(points, num_samples::Int)
    # Calculate distances between consecutive points
    distances = [norm(p2 .- p1) for (p1, p2) in zip(points[1:end-1], points[2:end])]
    total_distance = sum(distances)

    # Determine number of samples per segment proportionally
    segment_samples = [round(Int, d / total_distance * num_samples) for d in distances]

    # Ensure total samples match `num_samples`
    segment_samples[end] += num_samples - sum(segment_samples)

    # Perform linear interpolation for each segment
    interpolated_points = []
    input_indices = Int[]  # To store indices of input points in the output vector
    current_index = 1      # Tracks the position in the interpolated vector

    for i in 1:length(points)-1
        p1, p2 = points[i], points[i+1]
        n_samples = segment_samples[i]
        if i == 1  # Add the first input point only once
            push!(interpolated_points, p1)
            push!(input_indices, current_index)
            current_index += 1
        end
        # Linear interpolation for this segment
        for t in range(0.0, 1.0, length=n_samples + 1)[2:end-1]  # Avoid duplicating p1 or p2
            push!(interpolated_points, (1-t).*p1 .+ t.*p2)
            current_index += 1
        end
        # Add the second input point (end of the segment)
        push!(interpolated_points, p2)
        push!(input_indices, current_index)
        current_index += 1
    end

    return NTuple{length(points[1]),Float64}.(interpolated_points), input_indices
end



#Fourier Transforms
"""
    get_c_k(k::Tuple{Vararg{<:Real}}, c_iipDyn_mat::Array{T}, hte_lattice::Dyn_HTE_Lattice) -> T where {T}
    get_c_k(kvec::AbstractArray{<:Tuple{Vararg{<:Real}}}, c_iipDyn_mat::Array{T}, hte_lattice::Dyn_HTE_Lattice) -> Array{T} where {T}

Compute the spatial Fourier transform of dynamic correlators at momentum k.

This function calculates the Fourier transform of real-space correlators to reciprocal space,
assuming inversion symmetry to produce a real result. It sums over all basis states.

# Mathematical Background
The discrete Fourier transform is computed as:
```math
c(\\mathbf{k}) = \\frac{1}{N_{\\text{basis}}} \\sum_{b=1}^{N_{\\text{basis}}} \\sum_{i=1}^{N_{\\text{sites}}} \\cos(\\mathbf{k} \\cdot (\\mathbf{r}_i - \\mathbf{r}_{b})) \\, c_{ib}
```
where:
- ``\\mathbf{k}`` is the momentum vector
- ``\\mathbf{r}_i`` and ``\\mathbf{r}_b`` are the positions of site i and basis site b
- ``c_{ib}`` are the real-space correlators between sites i and basis site b
- ``N_{\\text{basis}}`` is the number of basis sites for normalization

# Arguments
- `k::Tuple{Vararg{<:Real}}` or `kvec::AbstractArray{<:Tuple{Vararg{<:Real}}}`: Single momentum vector or array of momentum vectors
- `c_iipDyn_mat::Array{T}`: Matrix of real-space dynamic correlators. See [`get_c_iipDyn_mat`](@ref).
- `hte_lattice::Dyn_HTE_Lattice`: Object containing lattice information. See [`Dyn_HTE_Lattice`](@ref).

# Returns
- Single k-point version: Value of the Fourier transform at momentum k with the same type as input correlators
- Multiple k-point version: Array of Fourier transforms with the same shape as `kvec`

# Examples
```julia
# Calculate correlator at a specific k-point K
k = (2π/3, 2π/sqrt(3))  # K-point in triangular lattice
c_kDyn = get_c_k(k, c_iipDyn_mat, hte_lattice)

# Create a path through high-symmetry points
path = [Γ, K, M, Γ]  # Common high-symmetry points in triangular lattice
k_vec, kticks_positions = create_brillouin_zone_path(path, 200)

# Calculate correlators along the path
c_k_vec = get_c_k(k_vec, c_iipDyn_mat, hte_lattice)

# Create a 2D grid of k-points for a BZ slice
N = 50
kx = (1:N)*2π/N
ky = (1:N)*2π/N
kmat = Tuple.([[1,1/sqrt(3)].*x .+ [1,-1/sqrt(3)].*y for x in kx, y in ky])
c_kDyn_grid = get_c_k(kmat, c_iipDyn_mat, hte_lattice)
```

# Notes
- Uses cosine transform to ensure real output assuming inversion symmetry
- Small values below 1e-12 are optionally set to zero (commented out in code)
- Result is normalized by the number of center sites
- Multiple k-point version uses broadcasting to efficiently apply the single k-point version
- Maintains the same shape/dimensionality as the input array when using the multiple k-point version

# See Also
- [`create_brillouin_zone_path`](@ref): Creates paths through high-symmetry points for band structure calculations
- [`inverse_fourier_transform`](@ref): Inverse transformation from k-space back to real space
- [`get_c_k_subl`](@ref): Sublattice-resolved version of this function
"""
function get_c_k(k::Tuple{Vararg{<:Real}},c_iipDyn_mat::Array{T},hte_lattice::Dyn_HTE_Lattice) where {T}

    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions

        if T == Taylor1{Float64}
            z = Taylor1(0)
        else
            z = zeros(size(c_iipDyn_mat[1]))
        end
        
        # Compute Fourier transformation at momentum k. The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
        for b in 1:length(lattice.unitcell.basis)
            for i in 1:length(lattice)
                z += cos(dot(k,getSitePosition(lattice,i).-getSitePosition(lattice,center_sites[b]))) *  c_iipDyn_mat[i,b]
            end
        end
        c_kDyn = z / length(center_sites) 
    
       #=  #set everything below 1e-12 to zero.
        for c_pos in  eachindex(c_kDyn)
            if abs(c_kDyn[c_pos]) < 1e-12
                c_kDyn[c_pos] = 0.0
            end
        end =#

    return return c_kDyn
end
function get_c_k(kvec::AbstractArray{<:Tuple{Vararg{<:Real}}},c_iipDyn_mat::Array{T},hte_lattice::Dyn_HTE_Lattice) where {T}
    fourier_transform(k) = get_c_k(k,c_iipDyn_mat,hte_lattice) 
    return fourier_transform.(kvec)
end

"""
    inverse_fourier_transform(kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}, 
                             c_kDyn::AbstractArray{T},
                             hte_lattice::Dyn_HTE_Lattice) -> Matrix{T} where {T}

Compute the inverse Fourier transform from reciprocal space to real space.

This function transforms k-space correlators back to real space, calculating correlations
between all lattice sites and basis sites.

# Mathematical Background

The inverse discrete Fourier transform is defined as:

```math
c_{i,b}= \\frac{1}{N_{\\mathbf{k}}} \\sum_{\\mathbf{k}} \\cos(\\mathbf{k} \\cdot (\\mathbf{r}_i - \\mathbf{r}_b)) \\, c(\\mathbf{k})
```

where:
- ``c_{i,b}`` is the real-space correlator between site ``i`` and basis site ``b``
- ``N_{\\mathbf{k}}`` is the total number of k-points for normalization
- ``\\mathbf{r}_i`` and ``\\mathbf{r}_b`` are the positions of site ``i`` and basis site ``b``
- ``c(\\mathbf{k})`` is the momentum-space correlator at wave vector ``\\mathbf{k}``

# Arguments
- `kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}`: Array of momentum vectors
- `c_kDyn::AbstractArray{T}`: Array of k-space correlators corresponding to `kvals`. Obtained from [`get_c_k`](@ref).
- `hte_lattice::Dyn_HTE_Lattice`: Object containing lattice information. See [`Dyn_HTE_Lattice`](@ref).

# Returns
- `Matrix{T}`: Matrix of real-space correlators with dimensions `[length(lattice), length(lattice.unitcell.basis)]`

# Examples
```julia
# Create a grid of k-points
N = 50
kx = (1:N)*2π/N
ky = (1:N)*2π/N
kmat = Tuple.([[1,1/sqrt(3)].*x .+ [1,-1/sqrt(3)].*y for x in kx, y in ky])

# Calculate k-space correlators
c_kDyn = get_c_k(kmat, c_iip_triang, triang_lattice)

# Transform back to real space
c_iip_triang == inverse_fourier_transform(kmat, c_kDyn, triang_lattice)
```

# Throws
- `ErrorException` if `c_kDyn` and `kvals` have different sizes

# Notes
- Uses cosine transform assuming inversion symmetry
- Result is normalized by the number of k-points

# See Also
- [`get_c_k`](@ref): Forward Fourier transform used to generate k-space input
- [`inverse_fourier_transform_subl`](@ref): Sublattice-resolved version of this function
"""
function inverse_fourier_transform(kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}
    ,c_kDyn::AbstractArray{T}
    ,
    hte_lattice::Dyn_HTE_Lattice)::Matrix{T} where {T}
    

    #check if kvals and c_kDyn have same dimensions
    if size(c_kDyn) != size(kvals)
        throw(error("c_kDyn and kvals have different sizes. They should be the same."))
    end

    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions
    (nx,ny) = size(kvals)

    c_iipDyn_mat = Array{T}(undef, length(lattice),length(lattice.unitcell.basis));
    #= Threads.@threads =# for k in 1:length(lattice)
        for b in 1:length(lattice.unitcell.basis)
        
            if T == Taylor1{Float64}
                z = Taylor1(0)
            elseif T == Float64
                z = 0.
            else
                z = zeros(size(c_kDyn[1]))
                end 

        for i in 1:nx,j in 1:ny
            # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
                z += cos(dot(kvals[i,j], getSitePosition(lattice,k).-getSitePosition(lattice,center_sites[b]))) *  c_kDyn[i,j]
        end
        c_iipDyn_mat[k,b] = z/length(eachindex(c_kDyn))
    end
    end
    return c_iipDyn_mat
end

#Sublattice Resolved Fourier Transforms
"""
    get_c_k_subl(k::Tuple{Vararg{<:Real}}, c_iipDyn_mat::Array{T}, hte_lattice::Dyn_HTE_Lattice) -> Matrix{Matrix{Float64}} where {T}
    get_c_k_subl(kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}, c_iipDyn_mat::Array{T}, hte_lattice::Dyn_HTE_Lattice) -> Array{Matrix{Matrix{Float64}}} where {T}

Compute the sublattice-resolved Fourier transform at momentum k.

This function calculates the sublattice-resolved Fourier transform, producing a matrix where
each element [b1,b2] represents correlations between sublattices b1 and b2.

# Mathematical Background

The sublattice-resolved Fourier transform is defined as:

```math
c_{b_1,b_2}(\\mathbf{k}) = \\sum_{i \\in b_2} e^{-i\\mathbf{k} \\cdot (\\mathbf{r}_i - \\mathbf{r}_{b_1})} \\, c_{i,b_1}
```

where:
- ``c_{b_1,b_2}(\\mathbf{k})`` is the sublattice-resolved correlator in momentum space
- ``b_1, b_2`` are sublattice indices  
- ``\\mathbf{r}_i`` is the position of site ``i`` belonging to sublattice ``b_2``
- ``\\mathbf{r}_{b_1}`` is the position of the center site of sublattice ``b_1``
- ``c_{i,b_1}`` is the real-space correlator between site ``i`` and the center site of sublattice ``b_1``

# Arguments
- `k::Tuple{Vararg{<:Real}}` or `kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}`: Single momentum vector or array of momentum vectors
- `c_iipDyn_mat::Array{T}`: Matrix of real-space dynamic correlators. See [`get_c_iipDyn_mat`](@ref).
- `hte_lattice::Dyn_HTE_Lattice`: Object containing lattice information. See [`Dyn_HTE_Lattice`](@ref).

# Returns
- Single k-point version: `Matrix{Matrix{Float64}}` with dimensions `[basis_size, basis_size]` where each element contains correlator coefficients
- Multiple k-point version: Array of sublattice-resolved Fourier transforms with the same shape as `kvals`

# Examples
```julia
# Calculate sublattice-resolved correlator for a single k-point
k_point = (2π/3, 2π/sqrt(3))  # K-point
c_k_subl = get_c_k_subl(k_point, c_iip_kagome, kagome_lattice)

# Shape of result is [basis_size, basis_size] matrices of coefficients
println("Sublattices: ", size(c_k_subl))  # For kagome, this is [3,3]

# Create a 2D grid of k-points
N = 30
kx = (1:N)*2π/N
ky = (1:N)*2π/N
kmat = Tuple.([[1,1/sqrt(3)].*x .+ [1,-1/sqrt(3)].*y for x in kx, y in ky])

# Calculate sublattice-resolved correlators for entire grid
c_kDyn = get_c_k_subl(kmat, c_iip_kagome, kagome_lattice)

invcorrs_kagome ≈ inverse_fourier_transform_subl(kmat, c_kDyn, kagome_lattice)
```

# Notes
- Uses complex exponential transform with e^(-i k·r) to allow for phase relationships between sublattices
- Takes the real part of the result assuming physical correlators should be real
- Uses `find_site_basis_label` to identify which sites belong to which sublattice
- Multiple k-point version uses multithreading with `@Threads.threads` for parallel computation
- Multiple k-point version returns a reshaped array matching the dimensions of the input `kvals`

# See Also
- [`get_c_k`](@ref): Non-sublattice-resolved version of this function
- [`inverse_fourier_transform_subl`](@ref): Inverse transformation from k-space back to real space
"""
function get_c_k_subl(k::Tuple{Vararg{<:Real}},c_iipDyn_mat::Array{T},hte_lattice::Dyn_HTE_Lattice) where {T}
    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions

    basis_size = length(lattice.unitcell.basis)
    label = find_site_basis_label(lattice)

    c_kDyn_mat = Array{Matrix{Float64}}(undef, basis_size,basis_size);
            
            # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
            for b1 in 1:basis_size
                for b2 in 1:basis_size
                
                    z = zeros(size(c_iipDyn_mat[1]))

                #find all indices that correspond to basis b2
                indexlist = findall(x->x==b2,label)

                # calculate correlator between the center site of b1 and all sites of b2 
                for index in indexlist
                    z += exp(-1im*dot(k, getSitePosition(lattice,index).-getSitePosition(lattice,center_sites[b1]))) *  c_iipDyn_mat[index,b1]
                end
                c_kDyn_mat[b1,b2] = real.(z)
                end

             end

    return c_kDyn_mat
end

function get_c_k_subl(kvals::AbstractArray{<:Tuple{Vararg{<:Real}}},c_iipDyn_mat::Array{T},hte_lattice::Dyn_HTE_Lattice) where {T}
    results = Array{Matrix{Matrix{Float64}}}(undef, length(kvals))

    @Threads.threads for i in eachindex(kvals)
        results[i] = get_c_k_subl(kvals[i], c_iipDyn_mat, hte_lattice)
    end
    return reshape(results, size(kvals)...)
end


"""
    inverse_fourier_transform_subl(kvals::AbstractArray{<:Tuple{Vararg{<:Real}}},
                                  c_kDyn_subl::Union{AbstractVector{Matrix{T}}, AbstractMatrix{Matrix{T}}, AbstractArray{Matrix{T}}},
                                  hte_lattice::Dyn_HTE_Lattice) -> Matrix{T} where {T}

Compute the inverse sublattice-resolved Fourier transform.

This function transforms sublattice-resolved k-space correlators back to real space,
calculating correlations between all lattice sites and basis sites.

# Mathematical Background

The inverse sublattice-resolved Fourier transform is defined as:

```math
c_{i,b_1} = \\frac{1}{N_{\\mathbf{k}}} \\sum_{\\mathbf{k}} e^{i\\mathbf{k} \\cdot (\\mathbf{r}_i - \\mathbf{r}_{b_1})} \\, c_{b_1,b_2}(\\mathbf{k})
```

where:
- ``c_{i,b_1}`` is the real-space correlator between site ``i`` and the center site of sublattice ``b_1``
- ``N_{\\mathbf{k}}`` is the total number of k-points for normalization
- ``\\mathbf{r}_i`` is the position of site ``i`` belonging to sublattice ``b_2``
- ``\\mathbf{r}_{b_1}`` is the position of the center site of sublattice ``b_1``
- ``c_{b_1,b_2}(\\mathbf{k})`` is the sublattice-resolved correlator in momentum space
- The sum is performed over all k-points and the result is filtered by sublattice membership

# Arguments
- `kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}`: Array of momentum vectors
- `c_kDyn_subl`: Array of sublattice-resolved k-space correlators corresponding to `kvals`
- `hte_lattice::Dyn_HTE_Lattice`: Object containing lattice information. See [`Dyn_HTE_Lattice`](@ref).

# Returns
- `Matrix{T}`: Matrix of real-space correlators with dimensions `[length(lattice), length(lattice.unitcell.basis)]`

# Examples
```julia
# Create a grid of k-points and calculate sublattice-resolved correlators
N = 30
kx = (1:N)*2π/N
ky = (1:N)*2π/N
kmat = Tuple.([[1,1/sqrt(3)].*x .+ [1,-1/sqrt(3)].*y for x in kx, y in ky])
c_kDyn = get_c_k_subl(kmat, c_iip_kagome, kagome_lattice)

# Transform back to real space
c_iip_kagome ≈ inverse_fourier_transform_subl(kmat, c_kDyn, kagome_lattice)
```

# Throws
- `ErrorException`: If `c_kDyn_subl` and `kvals` have different sizes

# Notes
- Uses complex exponential transform with e^(i k·r) for the inverse transform
- Takes the real part of the result assuming physical correlators should be real
- Result is normalized by the number of k-points
- Uses parallelization with `@Threads.threads` for the inner loop

# See Also
- [`get_c_k_subl`](@ref): Forward sublattice-resolved Fourier transform
- [`inverse_fourier_transform`](@ref): Non-sublattice-resolved inverse transform
"""
function inverse_fourier_transform_subl(kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}
    ,c_kDyn_subl::Union{
    AbstractVector{Matrix{T}},
    AbstractMatrix{Matrix{T}},
    AbstractArray{Matrix{T}}
    },
    hte_lattice::Dyn_HTE_Lattice)::Matrix{T} where {T}
   

    #check if kvals and c_kDyn_subl have same dimensions
    if size(c_kDyn_subl) != size(kvals)
        throw(error("c_kDyn_subl and kvals have different sizes. They should be the same."))
    end

    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions
    basis_size = length(lattice.unitcell.basis)
    label = find_site_basis_label(lattice)

    c_iipDyn_mat = Array{T}(undef, length(lattice),length(lattice.unitcell.basis));


    for b1 in 1:basis_size
            for b2 in 1:basis_size
                    indexlist = findall(x->x==b2,label)


    @Threads.threads for index in indexlist
                    if T == Taylor1{ComplexF64}||T == Taylor1{Float64}
                        z = Taylor1(0)
                        else
                        z = zeros(size(c_kDyn_subl[1][1]))
                    end 
                   

                    for i in eachindex(c_kDyn_subl)
                        z += exp(1im*dot(kvals[i], getSitePosition(lattice,index).-getSitePosition(lattice,center_sites[b1]))) *  c_kDyn_subl[i][b1,b2]
                    end
                  
                    c_iipDyn_mat[index,b1] = real(z/length(eachindex(c_kDyn_subl)))
                end
            end

        end
    return c_iipDyn_mat
end


###### moments, continued fractions and dynamical spin structure factors
"""
    get_moments_from_c_kDyn(c_kDyn::Matrix{Float64}) -> Vector{Polynomial{Float64}}

Extract frequency moments from dynamic correlator coefficients.

Computes the moments m(0), m(2), m(4), ..., m(2r_max) from the coefficients of the
dynamic correlator in k-space, applying the necessary sign flips to follow the mathematical
definition of moments.

# Arguments
- `c_kDyn::Matrix{Float64}`: Matrix of dynamic correlator coefficients in k-space

# Returns
- `Vector{Polynomial{Float64}}`: Vector of moment polynomials [m₀, m₂, m₄, ...]

# Examples
```julia
# Calculate moments for a k-point in a triangular lattice
k = (2π/3, 2π/sqrt(3))  # K-point
c_kDyn = get_c_k(k, c_iipDyn_mat, hte_lattice)
m_vec = get_moments_from_c_kDyn(c_kDyn)
```

# Notes
- The maximum moment order is determined by the size of the input matrix
- This implementation handles arbitrary-order moments up to r_max = floor(size(c_kDyn)[1]/2)

# See Also
- [`get_c_k`](@ref): Used to obtain the k-space correlators input to this function
- [`flipEvenIndexEntries`](@ref): Function used internally to apply alternating signs
- [`fromMomentsToδ`](@ref): Converts the resulting moments to continued fraction parameters
"""
function get_moments_from_c_kDyn(c_kDyn::Matrix{Float64})
   
    ### generalizes the following which only works for max_order=12 (r_max=6)
    #m0 = Polynomial(+flipEvenIndexEntries(c_kDyn_mat[:,1]))    
    #m2 = Polynomial(+flipEvenIndexEntries(c_kDyn_mat[3:end,2]))
    #m4 = Polynomial(-flipEvenIndexEntries(c_kDyn_mat[5:end,3]))
    #m6 = Polynomial(+flipEvenIndexEntries(c_kDyn_mat[7:end,4]))
    #m8 = Polynomial(-flipEvenIndexEntries(c_kDyn_mat[9:end,5]))
    #m10= Polynomial(+flipEvenIndexEntries(c_kDyn_mat[11:end,6]))
    #m12= Polynomial(-flipEvenIndexEntries(c_kDyn_mat[13:end,7]))
    #m_vec = [m0,m2,m4,m6,m8,m10,m12]

    r_max = Int(floor(size(c_kDyn)[1]/2))

    m0 = Polynomial(+flipEvenIndexEntries(c_kDyn[:,1]))

    m_vec = vcat(m0, [Polynomial(-(-1)^r*flipEvenIndexEntries(c_kDyn[(2*r+1):end,r+1])) for r in 1:r_max])
    
    return m_vec
end

 """
    fromMomentsToδ(m_vec::Vector{Float64}) -> (Vector{Float64}, Vector{Float64})
    fromMomentsToδ(m_vec::Vector{Polynomial{Float64, :x}}) -> Vector{Polynomial{Float64}}

Convert moments [m0, m2, m4, ...] to continued fraction parameters [δ0, δ1, δ2, ...].

This function implements the mathematical transformation from frequency moments to the
δ parameters used in continued fraction expansions for spectral functions. The resulting
δ parameters can be used with [`contFrac`](@ref) to evaluate the continued fraction
representation of dynamic correlators.

# Arguments
- `m_vec`: Vector of moments, either as scalars or polynomials in temperature (x=J/T)

# Returns
- For scalar moments: Tuple of (δ vector, r vector of indices [0,1,2,...])
- For polynomial moments: Vector of δ polynomials

# See Also
- [`contFrac`](@ref): Uses the δ parameters for continued fraction evaluation
- [`get_moments_from_c_kDyn`](@ref): Function that generates the input moments
- [`extrapolate_δvec`](@ref): Extends δ parameters to higher orders
- [`get_extrapolation_params`](@ref): Gets linear extrapolation parameters for terminators
"""
function fromMomentsToδ(m_vec::Vector{Float64})
   
    @assert length(m_vec)<=9

    m_vec_pad = 1*m_vec
    while length(m_vec_pad)<9
        append!(m_vec_pad,0.0)
    end
    m0,m2,m4,m6,m8,m10,m12,m14,m16 = m_vec_pad

    δ0 = m0
    
    δ1 = (m2/m0)
    
    δ2 = m4/m2-m2/m0
    
    δ3 = -(m0*(-m4^2 + m2*m6))/(m2^3 - m0*m2*m4)
    
    δ4 = ((m2*(m4^3 + m0*m6^2 + m2^2*m8 - m4*(2*m2*m6 + m0*m8)))/((m2^2 - m0*m4)*(-m4^2 + m2*m6)))
    
    δ5 = -(((m2^2 - m0*m4)*(-m6^3 + 2*m4*m6*m8 - m2*m8^2 + (-m4^2 + m2*m6)*m10))
    /((m4^2 - m2*m6)*(m4^3 + m0*m6^2 + m2^2*m8 - m4*(2*m2*m6 + m0*m8))))
    
    δ6 = (((m4^2 - m2*m6)*(-m6^4 - m4^2*m8^2 + m0*m8^3 + 2*m2*m4*m8*m10 - m2^2*m10^2 + m0*m4*m10^2 + (m4^3 + (m2^2 - m0*m4)*m8)*m12 + m6^2*(3*m4*m8 + 2*m2*m10 + m0*m12) 
    - 2*m6*((m4^2 + m0*m8)*m10 + m2*(m8^2 + m4*m12))))
    /((m4^3 + m0*m6^2 + m2^2*m8 - m4*(2*m2*m6 + m0*m8))*(m6^3 - 2*m4*m6*m8 + m2*m8^2 + (m4^2 - m2*m6)*m10)))
    
    δ7 = -(((m4^3 + m0*m6^2 + m2^2*m8 - m4*(2*m2*m6 + m0*m8))*(-m8^4 - m6^2*m10^2 + m2*m10^3 + 2*m4*m6*m10*m12 - m4^2*m12^2 + m2*m6*m12^2 
    + (m6^3 + (m4^2 - m2*m6)*m10)*m14 + m8^2*(3*m6*m10 + 2*m4*m12 + m2*m14) - 2*m8*((m6^2 + m2*m10)*m12 + m4*(m10^2 + m6*m14))))
    /((m6^3 - 2*m4*m6*m8 + m2*m8^2 + (m4^2 - m2*m6)*m10)*(m6^4 + m4^2*m8^2 - m0*m8^3 - 2*m2*m4*m8*m10 + m2^2*m10^2 - m0*m4*m10^2 
    - (m4^3 + (m2^2 - m0*m4)*m8)*m12 - m6^2*(3*m4*m8 + 2*m2*m10 + m0*m12) + 2*m6*((m4^2 + m0*m8)*m10 + m2*(m8^2 + m4*m12)))))
    
    δ8 = (((m6^3 - 2*m4*m6*m8 + m2*m8^2 + (m4^2 - m2*m6)*m10)*(m8^5 + m0*m10^4 - 2*m6^3*m10*m12 + 2*m2*m6*m10^2*m12 + 2*m0*m6*m10*m12^2 
    + m2^2*m12^3 - 2*m0*m6*m10^2*m14 + 2*m2*m6^2*m12*m14 - 2*m2^2*m10*m12*m14 + m0*m6^2*m14^2 + ((m6^2 - m2*m10)^2 - m0*m6^2*m12)*m16 - m8^3*(4*m6*m10 + 3*m4*m12 + 2*m2*m14 + m0*m16) 
    + m8^2*(3*m4*m10^2 + 3*m6^2*m12 + 4*m2*m10*m12 + m0*m12^2 + 4*m4*m6*m14 + 2*m0*m10*m14 + (m4^2 + 2*m2*m6)*m16) + m4^2*(m12*(m10^2 - 2*m6*m14) + 2*m6*m10*m16) 
    + m4^3*(m14^2 - m12*m16) + m8*(-3*m0*m10^2*m12 + 2*m4^2*m12^2 - 2*m6^3*m14 - 4*m4^2*m10*m14 - m0*m4*m14^2 + m0*m4*m12*m16 + 3*m6^2*(m10^2 - m4*m16) 
    + 2*m6*(-(m12*(m4*m10 + 2*m2*m12 + m0*m14)) + m0*m10*m16) - 2*m2*(m10^3 - m4*m12*m14 + m4*m10*m16) + m2^2*(m14^2 - m12*m16)) + m4*(2*m2*m10*(-m12^2 + m10*m14) 
    + m6^2*(m12^2 + 2*m10*m14) - m0*(m12^3 - 2*m10*m12*m14 + m10^2*m16) - 2*m6*(m10^3 + m2*(m14^2 - m12*m16)))))
    /((m6^4 + m4^2*m8^2 - m0*m8^3 - 2*m2*m4*m8*m10 + m2^2*m10^2 - m0*m4*m10^2 - (m4^3 + (m2^2 - m0*m4)*m8)*m12 - m6^2*(3*m4*m8 + 2*m2*m10 + m0*m12) 
    + 2*m6*((m4^2 + m0*m8)*m10 + m2*(m8^2 + m4*m12)))*(-m8^4 - m6^2*m10^2 + m2*m10^3 + 2*m4*m6*m10*m12 - m4^2*m12^2 + m2*m6*m12^2 
    + (m6^3 + (m4^2 - m2*m6)*m10)*m14 + m8^2*(3*m6*m10 + 2*m4*m12 + m2*m14) - 2*m8*((m6^2 + m2*m10)*m12 + m4*(m10^2 + m6*m14)))))
    
    δ_vec = [δ0,δ1,δ2,δ3,δ4,δ5,δ6,δ7,δ8][1:length(m_vec)] 
    return δ_vec , 1.0*collect(0:length(δ_vec)-1)
end
function fromMomentsToδ(m_vec::Vector{Polynomial{Float64, :x}})
    @assert length(m_vec)<=7

    m_vec_pad = 1*m_vec
    while length(m_vec_pad)<7
        append!(m_vec_pad,0.0)
    end
    m0,m2,m4,m6,m8,m10,m12 =m_vec_pad

    δ0 = m0
    
    δ1 = (m2//m0)
    
    δ2 = m4//m2-m2//m0
    
    δ3 = -(m0*(-m4^2 + m2*m6))//(m2^3 - m0*m2*m4)
    
    δ4 = ((m2*(m4^3 + m0*m6^2 + m2^2*m8 - m4*(2*m2*m6 + m0*m8)))//((m2^2 - m0*m4)*(-m4^2 + m2*m6)))
    
    δ5 = -(((m2^2 - m0*m4)*(-m6^3 + 2*m4*m6*m8 - m2*m8^2 + (-m4^2 + m2*m6)*m10))
    //((m4^2 - m2*m6)*(m4^3 + m0*m6^2 + m2^2*m8 - m4*(2*m2*m6 + m0*m8))))
    
    δ6 = (((m4^2 - m2*m6)*(-m6^4 - m4^2*m8^2 + m0*m8^3 + 2*m2*m4*m8*m10 - m2^2*m10^2 + m0*m4*m10^2 + (m4^3 + (m2^2 - m0*m4)*m8)*m12 + m6^2*(3*m4*m8 + 2*m2*m10 + m0*m12) 
    - 2*m6*((m4^2 + m0*m8)*m10 + m2*(m8^2 + m4*m12))))
    //((m4^3 + m0*m6^2 + m2^2*m8 - m4*(2*m2*m6 + m0*m8))*(m6^3 - 2*m4*m6*m8 + m2*m8^2 + (m4^2 - m2*m6)*m10)))
    

    δ_vec = [δ0,δ1,δ2,δ3,δ4,δ5,δ6][1:length(m_vec)]
    return δ_vec
end

"""
    contFrac(s::Number, δ_vec::Vector{Float64},a::Float64,b::Float64) -> Number

Evaluate an infinite continued fraction at complex frequency s using δ parameters and terminator function controlled by a,b.

This function recursively evaluates the continued fraction representation of the dynamic
correlator using the δ parameters from the moment expansion:

```math
CF(s) = \\cfrac{\\delta_0}{s + \\cfrac{\\delta_1}{s + \\cfrac{\\delta_2}{s + \\cfrac{\\delta_3}{s + \\ddots}}}}
```

For a finite number of δ parameters [δ₀, δ₁, ..., δₙ], the expansion becomes:
```math
CF(s) = \\cfrac{\\delta_0}{s + \\cfrac{\\delta_1}{s + \\cfrac{\\delta_2}{s + \\cfrac{\\ddots}{s + \\delta_n \\cdot T(s)}}}}
```
where T(s) is the terminator function provided by [`ContFracTerminator`](@ref).

It uses a terminator function that formally corresponds to a linear growth of the δ parameters beyond the last specified δ.
The linear growth is given by the equation δr = a*(r-(length(δ_vec)-1)) + b

# Arguments
- `s::Number`: Complex frequency at which to evaluate (typically s = iω + η)
- `δ_vec::Vector{Float64}`: Vector of δ parameters [δ₀, δ₁, ..., δᵣ]
- `a::Float64`: Slope of the linear growth for δ beyond the last specified δ
- `b::Float64`: Intercept of the linear growth for δ beyond the last specified δ

# See Also
- [`ContFracTerminator`](@ref): Provides the analytical terminator function T(s)
- [`fromMomentsToδ`](@ref): Converts moments to δ parameters used as input
- [`get_extrapolation_params`](@ref): Determines the terminator parameters a and b
- [`JSwithTerminator`](@ref): Uses continued fractions for dynamic structure factor computation
"""
function contFrac(s::Number,δ_vec::Vector{Float64},a::Float64,b::Float64)::Number
    
    if length(δ_vec)==1
        return  δ_vec[1]*ContFracTerminator(s,a,b)
    else
        return δ_vec[1]/(s+contFrac(s,δ_vec[2:end],a,b))
    end
end


"""
    extrapolate_δvec(δ_vec::Vector{Float64}, r_min::Int, r_max::Int, r_ext::Int, intercept0::Bool) -> Vector{Float64}

Extrapolate continued fraction parameters δ_vec=[δ₀, δ₁, ..., δᵣ] to higher orders using linear fitting.

This function fits a linear model to δ-parameters in the range r_min to r_max, then extrapolates 
additional parameters up to r_ext. The function supports fitting through the origin if needed.

# Arguments
- `δ_vec::Vector{Float64}`: Vector of existing δ parameters to extend
- `r_min::Int`: Minimum r index to use for fitting (inclusive)
- `r_max::Int`: Maximum r index to use for fitting (inclusive)
- `r_ext::Int`: Target maximum r index after extrapolation
- `intercept0::Bool`: If true, force fit line through origin; otherwise include intercept term

# Returns
- `Vector{Float64}`: Extended vector of δ parameters with length r_ext+1

# Examples
```julia
given a vector of moments m0_vec = [m₀, m₂, m₄, m₆, m₈]

# Convert moments to δ parameters and extrapolate
δ_vec, r_vec = fromMomentsToδ(m0_vec)
δ_vec_ext = extrapolate_δvec(δ_vec, r_max, r_max, r_ext, true)

```

# Notes
- Uses `LsqFit.curve_fit` to perform the linear regression
- If `intercept0=true`, uses model f(t) = a·t; otherwise uses f(t) = a·t + b
- The original δ parameters (δ₀ to δᵣₘₐₓ) are preserved in the output
- Typically used to extend δ parameters for continued fraction evaluation of dynamic correlators

# Throws
- `AssertionError` if r_max < r_min, r_ext ≤ r_max, or r_max+1 > length(δ_vec)

# See Also
- [`fromMomentsToδ`](@ref): Converts moments to δ parameters 
- [`get_extrapolation_params`](@ref): Alternative approach using parameters for terminator functions
- [`contFrac`](@ref): Uses extended δ parameters for continued fraction evaluation
"""
function extrapolate_δvec(δ_vec::Vector{Float64},r_min::Int,r_max::Int,r_ext::Int,intercept0::Bool)
    @assert r_max >= r_min
    @assert r_ext > r_max
    @assert r_max+1 <= length(δ_vec)
    ### define linear fit-function, fit and extrapolate
    if intercept0 
        function fa(t,p) return p[1] .* t end  
        p0 = [1.0]
        fit = LsqFit.curve_fit(fa, r_min:r_max, δ_vec[r_min+1:r_max+1], p0)
        return vcat(δ_vec[1:r_max+1],[fa(r,fit.param) for r in r_max+1:r_ext+1])
    else 
        function fab(t,p) return p[1] .* t .+ p[2] end
        p0 = [1.0,0.0]
        fit = LsqFit.curve_fit(fab, r_min:r_max, δ_vec[r_min+1:r_max+1], p0)
        return vcat(δ_vec[1:r_max+1],[fab(r,fit.param) for r in r_max+1:r_ext+1])
    end
end


"""
    get_extrapolation_params(δ_vec::Vector{Float64}, r_min::Int, r_max::Int, intercept0::Bool) -> Vector{Float64}

Get linear extrapolation parameters (a, b) for continued fraction terminator.

This function fits a linear model to δ-parameters in the range r_min to r_max and returns 
the parameters `a` (slope) and `b` (intercept) for use in [`ContFracTerminator`](@ref). These parameters
describe the asymptotic behavior δᵣ = a*(r - r_max) + b for r > r_max, which is then used
in [`contFrac`](@ref) to evaluate extrapolated continued fractions.

# Arguments
- `δ_vec::Vector{Float64}`: Vector of δ parameters [δ₀, δ₁, ..., δᵣ]
- `r_min::Int`: Minimum r index to use for fitting (inclusive)
- `r_max::Int`: Maximum r index to use for fitting (inclusive)
- `intercept0::Bool`: If true, force fit line through origin; otherwise include intercept term

# Returns
- `Vector{Float64}`: Two-element vector [a, b] where:
  - `a`: slope of the linear extrapolation
  - `b`: intercept parameter adjusted for continuation beyond r_max

# Notes
- Uses linear regression via `LsqFit.curve_fit` to determine the extrapolation model
- The returned parameters are specifically formatted for use with [`ContFracTerminator`](@ref)
- If `intercept0=true`, uses model f(t) = a·t; otherwise uses f(t) = a·t + b
- The `b` parameter is adjusted: b = intercept + (length(δ_vec)-1)*a

# Examples
```julia
# Get extrapolation parameters for continued fraction terminator
extrap_params = get_extrapolation_params(δ_vec, r_min, r_max, false)
a, b = extrap_params[1], extrap_params[2]

# Use in continued fraction evaluation
result = contFrac(s, δ_vec, a, b)
# or in Dynamic Structure Factor
JSwithTerminator(δ_vec[1:r_max+1], x, w, extrap_params)
```

# See Also
- [`ContFracTerminator`](@ref): Uses the returned parameters for analytical termination
- [`contFrac`](@ref): Main continued fraction evaluation function
- [`JSwithTerminator`](@ref): Dynamic structure factor computation using terminator
"""
function get_extrapolation_params(δ_vec::Vector{Float64},r_min::Int,r_max::Int,intercept0::Bool)
    @assert r_max >= r_min
    @assert r_max+1 <= length(δ_vec)
    if r_max==r_min && intercept0==false error("extrapolation of δ_vec not possible through just 1 point. Please set intercept0=true") end
    ### define linear fit-function, fit and extrapolate
    if intercept0 
        function fa(t,p) return p[1] .* t end  
        p0 = [1.0]
        fit = LsqFit.curve_fit(fa, r_min:r_max, δ_vec[r_min+1:r_max+1], p0)
        extrap_params =  [fit.param[1],0.0]  # return slope and intercept
    else 
        function fab(t,p) return p[1] .* t .+ p[2] end
        p0 = [1.0,0.0]
        fit = LsqFit.curve_fit(fab, r_min:r_max, δ_vec[r_min+1:r_max+1], p0)
        extrap_params = fit.param # return slope and intercept
    end
    a = extrap_params[1]
    b = extrap_params[2] + (length(δ_vec)-1)*a 

    return [a,b]
end



"""
    hermiteH(ν::Real, z::Number) -> Number
   
    #warning this function is only proportional to the Hermite function H_ν(z). It differs by a factor exp(zT^2)*(2^νT).

# See Also
- [`ContFracTerminator`](@ref): Uses this function for analytical continued fraction termination
"""
function hermiteH(ν::Real, z::Number)
    # promote to a common complex type
    T = Base.promote_type(typeof(ν + zero(ν)), typeof(z + zero(z)))
    νT = complex(T(ν))
    zT = complex(T(z))

    a1, b1 = -νT/2, T(1)/2
    a2, b2 = T(1)/2 - νT/2, T(3)/2

    #Factoring out exponential term
    M1 = #= exp(zT^2)* =#pFq((b1-a1,), (b1,), -zT^2)          # 1F1(-ν/2; 1/2; z^2)
    M2 = #= exp(zT^2)* =#pFq((b2-a2,), (b2,), -zT^2)          # 1F1(1/2 - ν/2; 3/2; z^2)

    return sqrt(pi)  *  #= (2^νT) *  =#( M1 / gamma((1 - ν)/2) - (2*zT) * M2 / gamma(-ν/2) )
end

"""
    ContFracTerminator(s::Number, a::Real, b::Real) -> Number

Exact expression for a continued fraction with linearly growing parameters.

This function provides the exact analytical solution for a continued fraction where the 
parameters follow δ(r) = a*r + b. The terminator is expressed in terms
of parabolic cylinder functions (Hermite functions) and provides a way to evaluate infinite
continued fractions with known asymptotic behavior.

# Arguments
- `s::Number`: Complex argument (typically frequency s = iω + η)
- `a::Real`: Slope parameter of the linear growth δ(r) = a*r + b
- `b::Real`: Intercept parameter of the linear growth

# Returns
- `Number`: Complex value representing the continued fraction terminator

# Notes
- This function is used internally by [`contFrac`](@ref)
- The parameters a and b are typically obtained from [`get_extrapolation_params`](@ref)
- Requires that a > 0 for convergence


# Mathematical Background
The infinite continued fraction with linearly growing parameters:
```math
\\Gamma_{a,b}(s) = \\cfrac{1}{s + \\cfrac{(1a + b)}{s + \\cfrac{(2a + b)}{s + \\cfrac{(3a + b)}{s + \\ddots}}}}
```
with a,b > 0 fulfills the functional equation:
```math
\\frac{1}{\\Gamma_{a,b}(s)} = s + (1a + b)\\Gamma_{a,1a+b}(s)
```

This is solved by a fraction of Hermite polynomials H_ν(z):
```math
\\Gamma_{a,b}(s) = \\frac{\\sqrt{2/a} \\, H_{-1-b/a}\\left(s/\\sqrt{2a}\\right)}{H_{-b/a}\\left(s/\\sqrt{2a}\\right)}
```

The special case b = 0 yields:
```math
\\Gamma_{a,0}(s) = \\sqrt{\\frac{\\pi}{2a}} e^{\\frac{s^2}{2a}} \\left[1 - \\text{erf}\\left(\\frac{s}{\\sqrt{2a}}\\right)\\right]
```
and analytically continues to a Gaussian:
```math
\\text{Re}[\\Gamma_{a,0}(i\\omega + 0^+)] = \\sqrt{\\frac{\\pi}{2a}} e^{-\\frac{\\omega^2}{2a}}
```


# See Also
- [`contFrac`](@ref): Uses this terminator for continued fraction evaluation
- [`get_extrapolation_params`](@ref): Determines the terminator parameters a and b
- [`JSwithTerminator`](@ref): Computes dynamic structure factors using this terminator
"""
function ContFracTerminator(s::Number,a::Real, b::Real)
     σ = a
     c = b
    return 1/2* sqrt(2/σ) * hermiteH(-(σ + c)/σ, s / sqrt(2*σ)) / hermiteH(-c/σ, s / sqrt(2*σ))
end

"""
    JSwithTerminator(δ_vec::Vector{Float64}, x::Float64, w::Float64, extrap_params::Vector{Float64}) -> Float64

Compute the dynamic structure factor J(S,ω) using continued fractions with terminator.

This function evaluates the dynamic structure factor (spectral function) at a given frequency
using the continued fraction representation with an analytical terminator. The terminator
parameters are obtained from linear extrapolation of the δ parameters, allowing for
evaluation of the infinite continued fraction.

# Arguments
- `δ_vec::Vector{Float64}`: Vector of continued fraction parameters [δ₀, δ₁, ..., δᵣ]
- `x::Float64`: Temperature-related parameter (typically β = 1/T)
- `w::Float64`: Frequency ω at which to evaluate the structure factor
- `extrap_params::Vector{Float64}`: Extrapolation parameters [a, b] from [`get_extrapolation_params`](@ref)

# Returns
- `Float64`: Value of the dynamic structure factor J(S,ω) at the specified frequency

# Mathematical Details
The function computes:
```math
J(S,\\omega) = \\frac{x\\omega}{1 - e^{-x\\omega}} \\cdot \\frac{1}{\\pi} \\cdot \\text{Re}[CF(i\\omega)]
```

where CF(iω) is the continued fraction evaluated at s = iω using [`contFrac`](@ref).

For the special cases x=0 or ω=0, it returns the simpler expression:
```math
J(S,\\omega) = \\frac{1}{\\pi} \\cdot \\text{Re}[CF(i\\omega)]
```

# Notes
- The extrapolation slope `a = extrap_params[1]` must be non-negative for convergence
- If a negative slope is detected, the function prints a warning and uses the absolute value

# Examples
```julia
# Get extrapolation parameters
extrap_params = get_extrapolation_params(δ_vec, r_min, r_max, false)

# Evaluate dynamic structure factor
x = 1.0  # inverse temperature
w = 2.0  # frequency
JS_value = JSwithTerminator(δ_vec, x, w, extrap_params)
```

# See Also
- [`get_extrapolation_params`](@ref): Determines the terminator parameters a and b
- [`contFrac`](@ref): Evaluates the continued fraction CF(iω)
- [`ContFracTerminator`](@ref): Provides the analytical terminator used internally
"""
function JSwithTerminator(δ_vec::Vector{Float64},x::Float64,w::Float64,extrap_params::Vector{Float64})::Float64
    a = extrap_params[1]
    b = extrap_params[2]
    if a < 0.0
        println("Slope of δ extrapolation must be non-negative")
        #throw(ErrorException("Slope of δ extrapolation must be non-negative"))
        a = abs(a)
    end

    res = 1/π * real(contFrac(1im * w ,δ_vec,a,b))
    if x==0.0 || w==0.0
        return res
    else
        return  x * w * 1/ (1 - exp(-x * w)) * res
    end
end 


"""
    get_JSkw_mat(method::String, x::Float64, k_vec::Vector, w_vec::Vector{Float64}, 
                c_iipDyn_mat::Array{Matrix{Rational{Int128}}}, lattice::Dyn_HTE_Lattice; 
                f::Float64=0.48, η::Float64=0.01, r_min::Int=3, r_max::Int=3, 
                r_ext::Int=1000, intercept0::Bool=true) -> Matrix{Float64}

Calculate dynamic spin structure factors J·S(k,ω) for given k-points and frequencies.

This function computes the dynamical structure factor matrix for a range of k-points
and frequencies. The bare series of the moments gets extrapolated with a specified method ("pade" or "u_pade") and the delta coeffficients 
of the continued fraction expansion get extrapolated linearly, according to the parameters "r_min", "r_max".


# Arguments
- `method::String`: Approximation method, either "pade" or "u_pade"
- `x::Float64`: Inverse temperature parameter x = J/T
- `k_vec::Vector`: Vector of k-points at which to calculate the structure factor
- `w_vec::Vector{Float64}`: Vector of frequency values ω/J
- `c_iipDyn_mat::Array{Matrix{Rational{Int128}}}`: Matrix of real-space dynamic correlators
- `lattice::Dyn_HTE_Lattice`: Lattice information object. See [`Dyn_HTE_Lattice`](@ref).
- `f::Float64=0.48`: Parameter for tanh transformation when using "u_pade" method
- `r_min::Int=3`: Minimum r value for δ parameter extrapolation
- `r_max::Int=3`: Maximum r value for δ parameter extrapolation
- `intercept0::Bool=true`: Whether to force linear extrapolation through origin

# Returns
- `Matrix{Float64}`: Matrix of J·S(k,ω) values with dimensions [length(k_vec), length(w_vec)]

# Examples
```julia
# define the lattice geometry
hte_lattice = getLattice(n_max,"chain")

#calculate the correlation function for all site combinations
c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs)

# define k and ω vectors 
k_vec = [(k*pi,0.0) for k in 0:0.1:2]
w_vec = collect(-3:0.1:3)

#calculate the spin structure factor for the given k and ω 
JSkw_mat = get_JSkw_mat("u_pade",4.0,k_vec,w_vec,c_iipDyn_mat,hte_lattice,r_min=2,r_max=3,f=0.48)


```

# Notes
- Issues warning if negative δ parameters are detected during extrapolation
- Progress is shown by printing k-point position during computation
- The "u_pade" method usually is more stable for small temperatures, but trivially gets constant for even smaller temperatures due to the saturation of tanh(f·x)

# See Also
- [`get_moments_from_c_kDyn`](@ref): Extracts frequency moments from dynamic correlator
- [`get_pade`](@ref): Constructs Padé approximants
- [`fromMomentsToδ`](@ref): Converts moments to δ parameters
- [`extrapolate_δvec`](@ref): Extrapolates δ parameters linearly
- [`get_c_k`](@ref): Fourier transform used internally
- [`get_LinearTrafoToCoeffs_u`](@ref): Matrix transformation for u-Padé method
"""
function get_JSkw_mat(method::String,x::Float64,k_vec::Vector,w_vec::Vector{Float64},c_iipDyn_mat::Array{Matrix{Rational{Int128}}},lattice::Dyn_HTE_Lattice;f::Float64=0.48,r_min::Int=3,r_max::Int=3,intercept0::Bool=true)

    JSkw_mat = 1.0*zeros(length(k_vec),length(w_vec))

    #pre-calculate the substitution matrix
    if method== "u_pade"
        substitution_matrix_arr = []
        for m_idx=1:6
            push!(substitution_matrix_arr, get_LinearTrafoToCoeffs_u(15-2*m_idx,f))
        end
    end


    for (k_pos,k) in enumerate(k_vec)
        println(k_pos,"/",length(k_vec))
        c_kDyn_mat = get_c_k([k],c_iipDyn_mat,lattice)[1]
        m_vec = get_moments_from_c_kDyn(c_kDyn_mat)[1:7]

        ###pade in x=J/T 
        if method=="pade"
            m_vec_extrapolated_pade = []
            for m_idx=1:length(m_vec)
                push!(m_vec_extrapolated_pade, get_pade(m_vec[m_idx],7-m_idx,7-m_idx))
            end
            δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec_extrapolated_pade])
        
        end



        ##pade with u=tanh(f*x) substitution
        if method == "u_pade"
            #if x= 0 we have to be careful with the substitution but case is trivial
            if x == 0
                m_vec_extrapolated_pade = []
                for m_idx=1:length(m_vec)
                    push!(m_vec_extrapolated_pade, get_pade(m_vec[m_idx],7-m_idx,7-m_idx))
                end
                δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec_extrapolated_pade])
            else

                m_vec_times_x =[m_vec[i]*Polynomial([0,1]) for i=1:length(m_vec)]
                m_vec_extrapolated_pade = []


                for m_idx=1:length(m_vec)-1
                    p_u = Polynomial(substitution_matrix_arr[m_idx]*coeffs(m_vec_times_x[m_idx]))
                    push!(m_vec_extrapolated_pade, get_pade(p_u,8-m_idx,7-m_idx))
                end
                
                δ_vec,r_vec = fromMomentsToδ([m(tanh(f*x))/x for m in m_vec_extrapolated_pade])
            end
        end


        ### exact extraplotation to r->infinity 

       
       # find last index where δ_vec is non-negative
        idx = findfirst(<(0), δ_vec)
        lastidx = isnothing(idx) ? length(δ_vec) : idx - 1
        r_max_eff = min(lastidx - 1, r_max)
        r_min_eff = min(r_min, r_max_eff)

        if r_max_eff < 1
            println("WARNING: negative δ1, putting δ0 = 0")
            extrap_params = [1.0,0.0]
            δ_vec = [0.0,1.0]
            r_max_eff = 0
        else
            extrap_params = get_extrapolation_params(δ_vec[1:r_max_eff+1],r_min_eff,r_max_eff,intercept0)
        end
        JSkw_mat[k_pos,:] = [JSwithTerminator(δ_vec[1:r_max_eff+1] ,x,w,extrap_params) for w in w_vec]
        
    end

    return JSkw_mat
end
