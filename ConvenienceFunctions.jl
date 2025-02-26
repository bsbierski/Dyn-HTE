using Symbolics, RobustPade, Polynomials, DifferentialEquations, LsqFit, TaylorSeries

###### get expansion coefficients for correlators G_ii' in real-space (dynamic-Matsubara and equal-time) 
function get_c_iipDyn_mat(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec)::Array{Matrix{Rational{Int64}}}
    """compute all non-trivial coefficients G_ii' on lattice L from the center_sites i to all other sites i' of the lattice"""
    GiipDyn_mat = Array{Matrix{Rational{Int64}}}(undef, lattice.length,length(lattice.unitcell.basis));
    Threads.@threads :dynamic for jp = 1:lattice.length
        println("bond "*string(jp)*" of "*string(lattice.length))
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

function get_c_iipEqualTime_mat(c_iipDyn_mat::Matrix{Matrix{Rational{Int64}}},max_order::Int)::Array{Rational{Int64}}
    """ perform frequency sum over real-space dynamic correlators to obtain equal time correlators """
    c_iipEqualTime_mat = Array{Rational{Int64}}(undef, length(c_iipDyn_mat[:,1]), length(c_iipDyn_mat[1,:]), max_order+1)
    for j in eachindex(c_iipDyn_mat[:,1])
        for b in eachindex(c_iipDyn_mat[1,:])
            c_iipEqualTime_mat[j,b,:] = [sum(c_iipDyn_mat[j,b][n+1,:] .* [1//1,1//12,1//720,1//30240,1//1209600,1//47900160,691//1307674368000,1//74724249600,3617//10670622842880000,43867//5109094217170944000]) for n in 0:max_order]
        end
    end
    return c_iipEqualTime_mat
end



###### bare series polynomial in Gii'(x,m) at Matsubara integer m truncated at n 
function get_TGiip_m_bare(c_iipDyn_mat::Matrix{Matrix{Rational{Int64}}},m::Int,n::Int)::Matrix{Polynomial}
    TGiip_bare = Array{Polynomial}(undef, lattice.length,length(lattice.unitcell.basis));
    for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
            if m==0
                TGiip_bare[jp,b] = Polynomial([1.0*c_iipDyn_mat[jp,b][nn+1,1]*(-1)^nn for nn in 0:n])
            else
                TGiip_bare[jp,b] = Polynomial([sum(c_iipDyn_mat[jp,b][nn+1,2:end] .* [1/(2*pi*m)^(2*l) for l in eachindex(c_iipDyn_mat[jp,b][1,2:end])])*(-1)^nn for nn in 0:n])
            end
        end
    end
    return TGiip_bare
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

function get_LinearTrafoToCoeffs_u(max_order::Int,f::Float64)::Matrix{Float64}
    """ get linear transform polynomial coeffs from x to u=tanh(fx) """ 
    """ to be used as res*coeffs_x """
    res = zeros(max_order+1,max_order+1)
    @variables x u
    x = taylor((atanh(u)/f), u, 0:max_order, rationalize=false)

    for n in 0:max_order
        xpn = simplify(x^n;expand=true)
        res[:,n+1] = Symbolics.value.(taylor_coeff(xpn,u,0:max_order,rationalize=false))
    end
    return res
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

    return interpolated_points, input_indices
end

function get_c_kDyn_mat(kvec,c_iipDyn_mat::Array{Matrix{Rational{Int64}}},lattice::Lattice,center_sites)::Vector{Matrix{Float64}}
    """computes the expansion coefficients c_k for the k-points in k_vec """
    BrillPath = Array{Matrix{Float64}}(undef,length(kvec));
    for i in eachindex(kvec)
        z = zeros(size(c_iipDyn_mat[1]))
        # Compute Fourier transformation at momentum k. The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
        for b in 1:length(lattice.unitcell.basis)
            for k in 1:length(lattice)
                z += cos(dot(kvec[i],getSitePosition(lattice,k).-getSitePosition(lattice,center_sites[b]))) *  c_iipDyn_mat[k,b]
            end
        end
       # z[z< 10^(-14)] .= 0
        BrillPath[i] = z 
    end
    
    return BrillPath
end



###### moments, continued fractions and dynamical spin structure factors
function get_moments_from_c_kDyn_mat(c_kDyn_mat::Matrix{Float64})
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

    r_max = Int(floor(size(c_kDyn_mat)[1]/2))

    m0 = Polynomial(+flipEvenIndexEntries(c_kDyn_mat[:,1]))

    m_vec = vcat(m0, [Polynomial(-(-1)^r*flipEvenIndexEntries(c_kDyn_mat[(2*r+1):end,r+1])) for r in 1:r_max])
    
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

function extrapolate_δvec(δ_vec::Vector{Float64},r_min::Int,r_max::Int,r_ext::Int,intercept0::Bool)
    """ extrapolate parameters of continued fraction δ_vec=[δ[0],δ[1],...,δ[R]] 
    using a linear interpolation for δ_vec[r_min] to δ[r_max], extrapolate δ[r_max+1]...δ[r_ext]. 
    If intercept0=true use line through origin. """
    @assert r_max > r_min
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

function get_JSkw_mat_x0(k_vec::Vector,w_vec::Vector{Float64},η::Float64,r_min::Int,r_max::Int,r_ext::Int,intercept0::Bool,c_iipDyn_mat::Array{Matrix{Rational{Int64}}},lattice::Lattice,center_sites)
    JSkw_mat = 1.0*zeros(length(k_vec),length(w_vec))

    for (k_pos,k) in enumerate(k_vec)
        c_kDyn_mat = get_c_kDyn_mat([k],c_iipDyn_mat,lattice,center_sites)[1]
        m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)

        ### pick x=0 for now (later check x>0, use PADE or IDA), get δs and plot them
        x=0.0
        δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec])

        δ_vec_ext = extrapolate_δvec(δ_vec,r_min,r_max,r_ext,intercept0)

        JSkw_mat[k_pos,:] = [JS(δ_vec_ext ,x,w,η) for w in w_vec]
    end

    return JSkw_mat
end


function get_JSkw_mat_finitex(diag_off_diag_flag,method::String,x::Float64,k_vec::Vector,w_vec::Vector{Float64},η::Float64,r_min::Int,r_max::Int,r_ext::Int,intercept0::Bool,c_iipDyn_mat::Array{Matrix{Rational{Int64}}},lattice::Lattice,center_sites)
    JSkw_mat = 1.0*zeros(length(k_vec),length(w_vec))
    max_order = Int((size(c_iipDyn_mat[1])[1]-1))

    for (k_pos,k) in enumerate(k_vec)
        #println(k_pos,"/",length(k_vec))
        c_kDyn_mat = get_c_kDyn_mat([k],c_iipDyn_mat,lattice,center_sites,diag_off_diag_flag)[1]
        m_vec = get_moments_from_c_kDyn_mat(c_kDyn_mat)[1:1+Int(floor(max_order/2))]

        ###PADE
        if method=="pade"
            m_vec_extrapolated_pade = Array{Any}(undef,length(m_vec))
            for m_idx=1:length(m_vec)
                m_vec_extrapolated_pade[m_idx] = get_pade(m_vec[m_idx],1+Int(floor(max_order/2))-m_idx,1+Int(floor(max_order/2))-m_idx)
            end
            δ_vec,r_vec = fromMomentsToδ([m(x) for m in m_vec_extrapolated_pade])
           println("Delta_vec (pade) for "*string(k)*"is:",δ_vec)
        end

        ###IDA
        if method =="ida"
            int_step_size = 0.01
            m_vec_extrapolated_ida = []
            for m_idx=1:2
                for idx=0:length(m_vec[m_idx])-1
                    if abs(m_vec[m_idx][idx])< 0.0000000001
                        m_vec[m_idx][idx] =0
                    end
                end
                IDA_parameters=[2,3-m_idx,3-m_idx] 
                IDA_approximant = get_intDiffApprox(m_vec[m_idx],collect(0:int_step_size:x+0.5),IDA_parameters[1],IDA_parameters[2],IDA_parameters[3])
                push!(m_vec_extrapolated_ida, Float64(IDA_approximant[round(Integer,1+1/int_step_size*x)]))
            end

            
            δ_vec,r_vec = fromMomentsToδ(Float64[float(m_vec_extrapolated_ida[i]) for i =1:length(m_vec_extrapolated_ida)])

            #println("Delta_vec (ida) for "*string(k)*"is:",δ_vec)
        end


        ###DIRECT DELTA EXTRAPOLATION VIA IDA
        if method=="directdelta"
            δ_vec_raw = fromMomentsToδ(m_vec)

            δ_vec = zeros(length(m_vec))

            for m_idx=1:2
                t = Taylor1(2*length(m_vec)-2)
                taylor_exp = δ_vec_raw[m_idx](t)
                taylor_coefficients = [taylor_exp[i] for i=0:length(taylor_exp)-1]

                taylor_poly = Polynomial(taylor_coefficients)
                
                #IDA extrapolation
                int_step_size = 0.01
                IDA_parameters=[2,3-m_idx,3-m_idx] 
                #println(taylor_poly)

                #println("m=",m_idx)
                IDA_approximant = get_intDiffApprox(taylor_poly,collect(0:int_step_size:x+0.5),IDA_parameters[1],IDA_parameters[2],IDA_parameters[3])

                δ_vec[m_idx] = IDA_approximant[round(Integer,1+1/int_step_size*x)]

                #pade extrapolation
                #δ_vec[m_idx] = get_pade(taylor_poly,6,6)(x)

            end

            #println("Delta_vec (direct delta) for "*string(k)*"is:",δ_vec)
        end


        ###DIRECT DELTA EXTRAPOLATION VIA IDA
        if method=="padedelta"
            δ_vec_raw = fromMomentsToδ(m_vec)

            δ_vec = zeros(length(m_vec))

            for δ_idx=1:length(δ_vec)
                t =  Taylor1(Float64, max_order+δ_idx)
                δ_poly= δ_vec_raw[δ_idx](t)
                δ_vec[δ_idx] = robustpade(δ_poly,6+δ_idx,6)(x)
            end

            println("Delta_vec (direct delta) for "*string(k)*"is:",δ_vec)
        end

        ###Now extrapolte deltas
        δ_vec_ext = extrapolate_δvec(δ_vec,r_min,r_max,r_ext,intercept0)

        JSkw_mat[k_pos,:] .= [JS(δ_vec_ext ,x,w,η) for w in w_vec]

        
    end

    return JSkw_mat
end


################################# BJÖRN STOPPED HERE ##################################

###### Brillouin Zone functions

function compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec)::Array{Matrix{Rational{Int64}}}
    """compute all correlations from the center_sites to all other sites of the lattice"""
    Correlators = Array{Matrix{Rational{Int64}}}(undef, lattice.length,length(lattice.unitcell.basis));
    Threads.@threads for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
        Correlators[jp,b] = mapreduce(permutedims, vcat, Calculate_Correlator_fast(LatGraph,center_sites[b],jp,max_order,gG_vec_unique,C_Dict_vec))
        end
    end
    return Correlators
end


function brillouin_zone_cut(kmat::Union{Matrix{Tuple{Float64,Float64}},Matrix{Tuple{Float64,Float64,Float64}}},Correlators::Matrix{Matrix{Rational{Int64}}},lattice::Lattice,center_sites)::Matrix{Matrix{Float64}}
    """computes the fourier transform along a 2D cut through the 2D or 3D k-space
        given the Correlation Matrix computet from compute_lattice_correlations """
    (nx,ny) = size(kmat)

    structurefactor = Array{Matrix{Float64}}(undef, nx,ny);
    Threads.@threads for i in 1:nx
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

function eval_correlator_LR_continuous_pad(Correlator,X,pade_order)
    """evaluate the correlator for imaginary frequencies by first fitting it to a continued fraction that preserves the continuity relation"""
    orderJ,orderω = size(Correlator)
    @variables x
    Gs = x-> sum(Correlator[i,1]*x^(i-1) for i =  1:orderJ) #static correlator
 
    Gs_pad = robustpade(Gs,pade_order[1],pade_order[2])

    G = Gs_pad(X)

    return G
   
end

function brillouin_zone_cut_Matrix(kmat::Union{Matrix{Tuple{Float64,Float64}},Matrix{Tuple{Float64,Float64,Float64}}},Correlators::Matrix{Matrix{Rational{Int64}}},lattice::Lattice,center_sites)::Matrix{Matrix{Matrix{Float64}}}
    """computes the fourier transform along a 2D cut through the 2D or 3D k-space
        given the Correlation Matrix computet from compute_lattice_correlations """
    (nx,ny) = size(kmat)

    structurefactor = Array{Matrix{Matrix{ComplexF64}}}(undef, nx,ny);
    Threads.@threads for i in 1:nx
        for j in 1:ny
            z = fourier_transform(kmat[i,j],Correlators,lattice,center_sites)
            structurefactor[i,j] = z #= / (length(lattice) * length(lattice.unitcell.basis)) =#
        end
    end
    return structurefactor
end

function fourier_transform(k,Correlators::Matrix{Matrix{Rational{Int64}}},lattice::Lattice,center_sites)::Matrix{Matrix{Float64}}
    """computes the fourier transform along a 2D cut through the 2D or 3D k-space
        given the Correlation Matrix computet from compute_lattice_correlations """
   
    basis_size = length(lattice.unitcell.basis)
    label = find_site_basis_label(lattice)


    structurefactor = Array{Matrix{ComplexF64}}(undef, basis_size,basis_size);
            
            # Compute Fourier transformation at momentum (kx, ky). The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
            for b1 in 1:basis_size
                for b2 in 1:basis_size
                
                    z = zeros(size(Correlators[1]))

                #find all indices that correspond to basis b2
                indexlist = findall(x->x==b2,label)

                # calculate correlator between the center site of b1 and all sites of b2 
                for index in indexlist
                    z += cos(dot(k, getSitePosition(lattice,index).-getSitePosition(lattice,center_sites[b1]))) *  Correlators[index,b1]
                end
                structurefactor[b1,b2] = z
                end

             end

    return structurefactor
end

function get_c_kDyn_mat(kvec,c_iipDyn_mat::Array{Matrix{Rational{Int64}}},lattice::Lattice,center_sites,diag_off_diag_flag::String)::Vector{Matrix{Float64}}
    """computes the expansion coefficients c_k for the k-points in k_vec """
    BrillPath = Array{Matrix{Float64}}(undef,length(kvec));

    label = find_site_basis_label(lattice)

    for i in eachindex(kvec)
        z = zeros(size(c_iipDyn_mat[1]))
        # Compute Fourier transformation at momentum k. The real-space position of the i-th spin is obtained via getSitePosition(lattice,i). 
        for b in 1:length(lattice.unitcell.basis)

            if diag_off_diag_flag == "total"
                indexlist = 1:length(lattice)
            elseif diag_off_diag_flag == "diag"
                indexlist = findall(x->x ==b,label)
            elseif diag_off_diag_flag == "off"
                indexlist = findall(x->x !=b,label)
   
            end

            for k in indexlist
                z += cos(dot(kvec[i],getSitePosition(lattice,k).-getSitePosition(lattice,center_sites[b]))) *  c_iipDyn_mat[k,b]
            end
        end
        z[abs.(z) .< 10^(-14)] .= 0
        BrillPath[i] = z 
    end
    return BrillPath
end
