using Symbolics, RobustPade, Polynomials, DifferentialEquations, LsqFit
using TaylorSeries



"""get expansion of real-space G_ii' between the sites in basis_positions and all sites in the Graph
    optional arguments:

    verbose: if true prints progress 
    max_order: restricts expansion to a maximum order of max_order.

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

"""get expansion of real-space G_ii' between the sites in hte_lattice.basis_positions and all sites in the hte_lattice.graph by using lattice symmetries
    optional arguments:

    verbose: if true prints progress 
    max_order: restricts expansion to a maximum order of max_order.

"""
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

function get_c_iipEqualTime_mat(c_iipDyn_mat::Matrix{Matrix{Rational{Int128}}})::Array{Vector{Rational{Int128}}}
    """ perform frequency sum over real-space dynamic correlators to obtain equal time correlators """
    max_order_plus1 = size(c_iipDyn_mat[1,1])[1]
    c_iipEqualTime_mat = Array{Vector{Rational{Int64}}}(undef, length(c_iipDyn_mat[:,1]), length(c_iipDyn_mat[1,:]))
    for j in eachindex(c_iipDyn_mat[:,1])
        for b in eachindex(c_iipDyn_mat[1,:])
            c_iipEqualTime_mat[j,b] = [sum(c_iipDyn_mat[j,b][n,:] .* [1//1,1//12,1//720,1//30240,1//1209600,1//47900160,691//1307674368000,1//74724249600,3617//10670622842880000,43867//5109094217170944000]) for n in 1:max_order_plus1]
        end
    end
    return c_iipEqualTime_mat
end



###### bare series polynomial in Gii'(x,m) at Matsubara integer m truncated at n 

""" get the expansion of the Matsubara correlator TGii'(iνm) as x-Polyomial for spatial entries i,ip of c_iipDyn_mat"""
function get_TGiip_Matsubara_xpoly(c_iipDyn_mat::Matrix{Matrix{Rational{Int128}}},i::Int,ip::Int,m::Int)
    if m==0
        p_x = 1.0*Polynomial(flipEvenIndexEntries(c_iipDyn_mat[i,ip][:,1]))
    else
        coeffs_m = [sum([c_iipDyn_mat[i,ip][n+1,lhalf+1] * 1/(2*π*m)^(2*lhalf) for lhalf in 1:9]) for n in 0:n_max]
        p_x = 1.0*Polynomial(flipEvenIndexEntries(coeffs_m))
    end

    return p_x
end


function flipEvenIndexEntries(v)
    """ v=[a,b,c,d,...] -> [+a,-b,+c,-d,...] """
    signs = [-1*(-1)^n for n in eachindex(v)]
    return v .* signs
end

###### resummation tools for polynomial p

function get_pade(p::Polynomial,N::Int,M::Int)
    """ Padé approximant (use RobustPade to avoid dividing by zero and just return 0//1) """
    return robustpade(p,N,M)
end

function get_intDiffApprox(p::Polynomial,x_vec::Vector{Float64},M::Int,L::Int,N::Int)
    """ Integrated differential approximant, setup of the ODE and return solution at x_vec """
    @assert M+L+N+2 <= Polynomials.degree(p)
    pp= Polynomials.derivative(p)
    @variables x

    f = Symbolics.series(p.coeffs,x)
    fp =Symbolics.series(pp.coeffs,x)

    cs, = @variables c[1:(M+L+N+2)]
    Q = Symbolics.series([c[k] for k in 1:M+1],x)
    P = expand(1.0 + x*Symbolics.series([c[k] for k in M+2:M+1+L],x))
    R = Symbolics.series([c[k] for k in M+2+L:M+L+N+2],x)
    s = Q * fp + P * f + R
    eqns = [taylor_coeff(s,x,m) ~ 0 for m in 0:M+L+N+1]
    res = symbolic_linear_solve(eqns, cs)

    ## solve differential equation
    function Q_fit(x) return Symbolics.series([res[k] for k in 1:M+1],x) end
    function P_fit(x) return expand(1.0 + x*Symbolics.series([res[k] for k in M+2:M+1+L],x)) end
    function R_fit(x) return Symbolics.series([res[k] for k in M+2+L:M+L+N+2],x) end

    g(f, ppp, x) = -(P_fit(x)*f+R_fit(x))/(Q_fit(x))
    f0 = p(0.0)
    xspan = (0.0, maximum(x_vec))
    prob = ODEProblem(g, f0, xspan)
    sol = DifferentialEquations.solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8, saveat=x_vec)
    return sol.u
end

###### variable transform from x to u=tanh(fx)
function get_p_u(coeffs_x::Vector{Float64},f::Float64)
    """ transform polynomial in x (defined via coeffs_x) to polynomial in u=tanh(fx) truncated to degree length(coeffs_x)-1 """ 
    """ this is too slow for transforming many coeffs_x, use the linear trafo instead """
    @variables x u
    x = taylor(atanh(u)/f, u, 0:(length(coeffs_x)-1), rationalize=false)
    p_u_ext = simplify(series(coeffs_x,x);expand=true)
    p_u = Polynomial(Symbolics.value.(taylor_coeff(p_u_ext,u,0:12,rationalize=false)),:u)
    return p_u
end

function get_LinearTrafoToCoeffs_u(max_order::Int, f::Float64)::Matrix{Float64}
""" get linear transform polynomial coeffs from x to u=tanh(fx) """ 
    """ to be used as res*coeffs_x """
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


###### k-space functions

function create_brillouin_zone_path(points, num_samples::Int)
    """ create a linear interpolation between an arbitrary number of (high symmetry) points in BZ """
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
function get_c_k(k::Tuple{Vararg{<:Real}},c_iipDyn_mat::Array{T},hte_lattice::Dyn_HTE_Lattice) where {T}
    """ computes the spatial FT of c_iipDyn for momentum k """
    """ assumes inversion symmetry of the lattice to get real FT transform """
    """ sums over all basis states:  """

    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions
        z = zeros(size(c_iipDyn_mat[1]))
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


function get_c_k(kvec::AbstractArray{<:Tuple{Vararg{<:Real}}}
    ,c_iipDyn_mat::Array{T},hte_lattice::Dyn_HTE_Lattice) where {T}

        fourier_transform(k) = get_c_k(k,c_iipDyn_mat,hte_lattice) 

        return fourier_transform.(kvec)
end


function inverse_fourier_transform(kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}
    ,c_kDyn::AbstractArray{T}
    ,
    hte_lattice::Dyn_HTE_Lattice)::Matrix{T} where {T}
    """computes the inverse fourier transform for sublattice resolved fourier transforms"""

    #check if kvals and c_kDyn have same dimensions
    if size(c_kDyn) != size(kvals)
        throw(error("c_kDyn and kvals have different sizes. They should be the same."))
    end

    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions
    (nx,ny) = size(kmat)

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
                z += cos(dot(kmat[i,j], getSitePosition(lattice,k).-getSitePosition(lattice,center_sites[b]))) *  c_kDyn[i,j]
        end
        c_iipDyn_mat[k,b] = z/length(eachindex(c_kDyn))
    end
    end
    return c_iipDyn_mat
end

#Sublattice Resolved Fourier Transforms

function get_c_k_subl(k::Tuple{Vararg{<:Real}},c_iipDyn_mat::Array{T},hte_lattice::Dyn_HTE_Lattice) where {T}
    """ computes the sublattice resolved spatial FT of c_iipDyn for momentum k """
    """ assumes inversion symmetry of the lattice to get real FT transform """
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

function get_c_k_subl(
    kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}
    ,c_iipDyn_mat::Array{T},hte_lattice::Dyn_HTE_Lattice) where {T}

    fourier_transform(k) = get_c_k_subl(k,c_iipDyn_mat,hte_lattice) 

    return fourier_transform.(kvals)
end

function inverse_fourier_transform_subl(kvals::AbstractArray{<:Tuple{Vararg{<:Real}}}
    ,c_kDyn_subl::Union{
    AbstractVector{Matrix{T}},
    AbstractMatrix{Matrix{T}},
    AbstractArray{Matrix{T}}
    },
    hte_lattice::Dyn_HTE_Lattice)::Matrix{T} where {T}
    """computes the inverse fourier transform for sublattice resolved fourier transforms"""

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


                for index in indexlist
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
function get_moments_from_c_kDyn(c_kDyn::Matrix{Float64})
    """ get the moments m(0),m(2)),m(4),...,m(2r_max) from the coefficients c_kDyn
    (need to flip the even indices in front of x^odd to comply with definition of c)v"""

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

function fromMomentsToδ(m_vec::Vector{Float64})
    """ convert [m0,m2,m4,...] of any length <9 to [δ0,δ1,...] of the same length, also return r_vec=[0,1,2,...] """
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


###Ruben? Kann das Weg?
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


function contFrac(s::Number,δ_vec::Vector{Float64})::Number
    """ continued fraction in variable s using δ_vec=[δ0,δ1,...,δr] and r-pole termination time τ"""
    if length(δ_vec)==1
        return  abs(δ_vec[1])^0.5
    else
        return δ_vec[1]/(s+contFrac(s,δ_vec[2:end]))
    end
end


#Björn überarbeiten
function extrapolate_δvec(δ_vec::Vector{Float64},r_min::Int,r_max::Int,r_ext::Int,intercept0::Bool)
    """ extrapolate parameters of continued fraction δ_vec=[δ[0],δ[1],...,δ[R]] 
    using a linear interpolation for δ_vec[r_min] to δ[r_max], extrapolate δ[r_max+1]...δ[r_ext]. 
    If intercept0=true use line through origin. """
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

function JS(δ_vec::Vector{Float64},x::Float64,w::Float64,η::Float64)::Float64
    """ get dynamical spin structure factor (times J) from δ_wec at w=ω/J and broadening η"""
    res = 1/π * real(contFrac(1im * w + η,δ_vec))
    if x==0.0 || w==0.0
        return res
    else
        return  x * w * 1/ (1 - exp(-x * w)) * res
    end
end 



function get_JSkw_mat(method::String,x::Float64,k_vec::Vector,w_vec::Vector{Float64},c_iipDyn_mat::Array{Matrix{Rational{Int128}}},lattice::Dyn_HTE_Lattice;f::Float64=0.48,η::Float64=0.01,r_min::Int=3,r_max::Int=3,r_ext::Int=1000,intercept0::Bool=false)
    """ get the dynamical spin structure factor from the correlation matrix c_iipDyn_mat 
    using pade approximants for the moments either in the variable x = J/T ("pade") or in the variable
    u = tanh(f*x) ("u_pade")   """


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


        ###Now extrapolte deltas
        δ_vec_ext = extrapolate_δvec(δ_vec,r_min,r_max,r_ext,intercept0)


        #if deltas have negative slope give warning
        if δ_vec_ext[end] <0
            δ_vec_ext = extrapolate_δvec(δ_vec,2,2,r_ext,false)
            println("WARNING: NEGATIVE δ parameters (pade might fail)")
            
        else
            JSkw_mat[k_pos,:] = [JS(δ_vec_ext ,x,w,η) for w in w_vec]
        end

        
    end

    return JSkw_mat
end



function extrapolate_series(series,method::String,parameters)
    
    if method == "pade"
        return get_pade(series,parameters[1],parameters[2])
    elseif method == "u_pade"
        substitution_matrix = get_LinearTrafoToCoeffs_u(length(coeffs(series))-1,parameters[3])
        p_u = Polynomial(substitution_matrix*coeffs(series))
        return get_pade(p_u,parameters[1],parameters[2])
    end

end



"""
Find the smallest x in [x_min, x_max] such that |f1(x) - f2(x)| > epsilon.
Returns the x value or `nothing` if no such point exists.
"""
function find_divergence_point(f1, f2, epsilon; x_min=0.0, x_max=10.0, step=0.01)
    diff(x) = abs(f1(x) - f2(x))

    x = x_min
    while x + step <= x_max
        if diff(x) ≤ epsilon && diff(x + step) > epsilon
            # Define a one-argument function for root finding
            return (x + step / 2)
        end
        x += step
    end

    return nothing  # No divergence found in the interval
end
