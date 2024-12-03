using Symbolics, RobustPade, Polynomials, DifferentialEquations

###### get expansion coefficients for correlators G_ii' in real-space (dynamic-Matsubara and equal-time) 
function get_c_iipDyn_mat(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec)::Array{Matrix{Rational{Int64}}}
    """compute all non-trivial coefficients G_ii' on lattice L from the center_sites i to all other sites i' of the lattice"""
    GiipDyn_mat = Array{Matrix{Rational{Int64}}}(undef, lattice.length,length(lattice.unitcell.basis));
    Threads.@threads for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
            GiipDyn_mat[jp,b] = mapreduce(permutedims, vcat, Calculate_Correlator_fast(LatGraph,center_sites[b],jp,max_order,gG_vec_unique,C_Dict_vec))
        end
    end
    return GiipDyn_mat
end
function get_c_iipDyn_mat_slow(LatGraph,lattice,center_sites,max_order,gG_vec,C_Dict_vec)::Array{Matrix{Rational{Int64}}}
    """slow version of get_c_iipDyn_mat not using gG_vec_unique but gG_vec (legacy, for chain)"""
    c_iipDyn_mat = Array{Matrix{Rational{Int64}}}(undef, lattice.length,length(lattice.unitcell.basis));
    for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
            c_iipDyn_mat[jp,b] = mapreduce(permutedims, vcat, Calculate_Correlator(LatGraph,center_sites[b],jp,max_order,gG_vec,C_Dict_vec))
        end
    end
    return c_iipDyn_mat
end

function get_c_iipEqualTime_mat(GiipDyn_mat::Matrix{Matrix{Rational{Int64}}},max_order::Int)::Matrix{Rational{Int64}}
    """ perform frequency sum over real-space dynamic correlators to obtain equal time correlators """
    GiipEqualTime_mat = Matrix{Rational{Int64}}(undef, length(GiipDyn_mat),max_order+1)
    for j in eachindex(GiipDyn_mat)
        GiipEqualTime_mat[j,:] = [sum(GiipDyn_mat[j][n+1,:] .* [1//1,1//12,1//720,1//30240,1//1209600,1//47900160,691//1307674368000,1//74724249600,3617//10670622842880000,43867//5109094217170944000]) for n in 0:max_order]
    end
    return GiipEqualTime_mat
end

###### bare series polynomial in Gii'(x,m) at Matsubara integer m truncated at n 
function get_TGiip_m_bare(c_iipDyn_mat::Matrix{Matrix{Rational{Int64}}},m::Int,n::Int)::Matrix{Polynomial}
    TGiip_bare = Array{Polynomial}(undef, lattice.length,length(lattice.unitcell.basis));
    for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
            if m==0
                TGiip_bare[jp,b] = Polynomial([1.0*c_iipDyn_mat[jp,b][nn+1,1] for nn in 0:n])
            else
                TGiip_bare[jp,b] = Polynomial([sum(c_iipDyn_mat[jp,b][nn+1,2:end] .* [1/(2*pi*m)^(2*l) for l in eachindex(c_iipDyn_mat[jp,b][1,2:end])]) for nn in 0:n])
            end
        end
    end
    return TGiip_bare
end

###### resummation tools
function get_pade(p::Polynomial,N::Int,M::Int)
    """ Padé approximant (use RobustPade to avoid dividing by zero and just return 0//1) """
    return robustpade(p,N,M)
end

function get_intDiffApprox(p::Polynomial,x_vec::Vector{Float64},M::Int,L::Int,N::Int)
    """ Integrated differential approximant, setup of the ODE and return solution at x_vec """
    @assert M+L+N+2 <= Polynomials.degree(p)
    pp= derivative(p)
    @variables x

    f=series(p.coeffs,x)
    fp=series(pp.coeffs,x)

    cs, = @variables c[1:(M+L+N+2)]
    Q = series([c[k] for k in 1:M+1],x)
    P = expand(1.0 + x*series([c[k] for k in M+2:M+1+L],x))
    R = series([c[k] for k in M+2+L:M+L+N+2],x)
    s = Q * fp + P * f + R
    eqns = [taylor_coeff(s,x,m) ~ 0 for m in 0:M+L+N+1]
    res = symbolic_linear_solve(eqns, cs)

    ## solve differential equation
    function Q_fit(x) return series([res[k] for k in 1:M+1],x) end
    function P_fit(x) return expand(1.0 + x*series([res[k] for k in M+2:M+1+L],x)) end
    function R_fit(x) return series([res[k] for k in M+2+L:M+L+N+2],x) end

    g(f, ppp, x) = -(P_fit(x)*f+R_fit(x))/(Q_fit(x))
    f0 = p(0.0)
    xspan = (0.0, maximum(x_vec))
    prob = ODEProblem(g, f0, xspan)
    sol = DifferentialEquations.solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8, saveat=x_vec)
    return sol.u
end






################################# BJÖRN STOPPED HERE ##################################


###### Brillouin Zone functions
function brillouin_zone_cut(kmat::Union{Matrix{Tuple{Float64,Float64}},Matrix{Tuple{Float64,Float64,Float64}}},Correlators::Matrix{Matrix{Rational{Int64}}},lattice::Lattice,center_sites)::Matrix{Matrix{Float64}}
    """computes the fourier transform along a 2D cut through the 2D or 3D k-space
        given the Correlation Matrix computet from compute_lattice_correlations """
    (nx,ny) = size(kmat)

    structurefactor = Array{Matrix{Float64}}(undef, nx,ny);
    for i in 1:nx
    for j in 1:ny
            z = zeros(size(Correlators[1]))
            # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
            for b in 1:length(lattice.unitcell.basis)
                for k in 1:length(lattice)
                    z += cos(dot(kmat[i,j], getSitePosition(lattice,k).-getSitePosition(lattice,center_sites[b]))) *  Correlators[k,b]
                end
            end
            structurefactor[i,j] = z #= / (length(lattice) * length(lattice.unitcell.basis)) =#
        end
    end
    return structurefactor
end

function brillouin_zone_path(kvec,Correlators::Matrix{Matrix{Rational{Int64}}},lattice::Lattice,center_sites)::Vector{Matrix{Float64}}
    """computes the fourier transform along a 1D path through k-space given the path computed by create_brillouin_zone_path """
    BrillPath = Array{Matrix{Float64}}(undef,length(kvec));
    for i in eachindex(kvec)
        z = zeros(size(Correlators[1]))
        # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
        for b in 1:length(lattice.unitcell.basis)
            for k in 1:length(lattice)
                z += cos(dot(kvec[i],getSitePosition(lattice,k).-getSitePosition(lattice,center_sites[b]))) *  Correlators[k,b]
            end
        end
        BrillPath[i] = z 
    end
    return BrillPath
end

function create_brillouin_zone_path(points, num_samples::Int)
    """Creates a linear interpolation between high symmetry points"""
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

    return interpolated_points, input_indices
end

#### evaluation functions
function eval_correlator_LR_continuous_pad(Correlator,ω,JoverT,pade_order)
    """evaluate the correlator for real frequencies by first fitting it to a continued fraction that preserves the continuity relation """
    orderJ,orderω = size(Correlator)

    @variables x
    δ = 0.1
    cs = [ x-> sum([Correlator[i,n]*x^(i-1) for i =  1:orderJ#= (1+max_order) =#]) for n = 2:orderω]
    Gs = x-> sum(Correlator[i,1]*x^(i-1) for i =  1:orderJ) #static correlator
 
    #pade aprpox

    cs_pad = [robustpade(c,8,0) for c in cs]
    Gs_pad = robustpade(Gs,pade_order[1],pade_order[2])

    c = [cf(JoverT) for cf in cs_pad]
    G = Gs_pad(JoverT)
    return  c[1] / (-(JoverT*(ω + 1im*δ))^2 - c[2]/c[1] - (-c[2]^2 + c[1]*c[3]) /
            (c[1]^2 * (-(JoverT*(ω + 1im*δ))^2 + (G * (c[2]^2 - c[1]*c[3])) /
            (c[1]^3 + G*c[1]*c[2]))))
end



#= 
### this function is unstable
 function eval_correlator_LR_continuous_pad(Correlator,ω,JoverT,pade_order)
    orderJ,orderω = size(Correlator)

    @variables x
    δ = 0.3
    cs = [ x-> sum([Correlator[i,n]*x^(i-1) for i =  1:orderJ#= (1+max_order) =#]) for n = 2:orderω]
    Gs = x-> sum(Correlator[i,1]*x^(i-1) for i =  1:orderJ) #static correlator

    #pade aprpox

    cs_pad = [robustpade(c,pade_order[1],pade_order[2]) for c in cs]
    Gs_pad = robustpade(Gs,pade_order[1],pade_order[2])

    c = [cf(JoverT) for cf in cs]
    G = Gs_pad(JoverT)
    numerator = c[1]
    denominator = (1im*JoverT*(ω + 1im*δ))^2 - c[2]/c[1] - (-c[2]^2 + c[1]*c[3]) / 
    (c[1]^2 * ((1im*JoverT*(ω + 1im*δ))^2 + (-c[2]^3 + 2*c[1]*c[2]*c[3] - c[1]^2*c[4]) / 
    (c[1]*(-c[2]^2 + c[1]*c[3])) - 
    (c[1] * (-c[3]^3 + 2*c[2]*c[3]*c[4] - c[1]*c[4]^2 - c[2]^2*c[5] + 
    c[1]*c[3]*c[5])) / 
    ((c[2]^2 - c[1]*c[3])^2 * ((1im*JoverT*(ω + 1im*δ))^2 + 
    (c[1]*(c[1] + (G*c[2])/c[1]) * 
    (-c[3]^3 + 2*c[2]*c[3]*c[4] - c[1]*c[4]^2 - c[2]^2*c[5] + c[1]*c[3]*c[5])) / 
    ((c[2]^2 - c[1]*c[3])^2 * 
    ((G*(-c[2]^2 + c[1]*c[3]))/c[1]^2 + 
    (-c[2]^3 + 2*c[1]*c[2]*c[3] - c[1]^2*c[4])/(-c[2]^2 + c[1]*c[3]) + 
    (G*c[2]*(-c[2]^3 + 2*c[1]*c[2]*c[3] - c[1]^2*c[4])) / 
    (c[1]^2*(-c[2]^2 + c[1]*c[3]))))))))

    return numerator / denominator

end
=#
 
 
function eval_correlator_LR_continuous_pad_Mats(Correlator,ω,X,pade_order)
    """evaluate the correlator for imaginary frequencies by first fitting it to a continued fraction that preserves the continuity relation"""
    orderJ,orderω = size(Correlator)

    @variables x
    cs = [ x-> sum([Correlator[i,n]*x^(i-1) for i =  1:orderJ#= (1+max_order) =#]) for n = 2:orderω]
    Gs = x-> sum(Correlator[i,1]*x^(i-1) for i =  1:orderJ) #static correlator
 
    #pade aprpox

    cs_pad = [robustpade(c,orderJ-1,0) for c in cs]
    Gs_pad = robustpade(Gs,pade_order[1],pade_order[2])

    c = [cf(X) for cf in cs_pad]
    G = Gs_pad(X)

    if ω ==0 
        return G
    end
    
    return  c[1] / (-(1im*ω)^2 - c[2]/c[1] - (-c[2]^2 + c[1]*c[3]) /
            (c[1]^2 * (-(1im*ω)^2 + (G * (c[2]^2 - c[1]*c[3])) /
            (c[1]^3 + G*c[1]*c[2]))))
end
 


