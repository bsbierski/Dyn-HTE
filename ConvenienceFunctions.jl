#using RobustPade
using Symbolics

###### get Dyn-HTE correlators in real-space 
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

function getCorrelators_equalTime(Correlators::Matrix{Matrix{Rational{Int64}}},max_order::Int)::Matrix{Rational{Int64}}
    """ perform frequency sum over real-space dynamic correlators to obtain equal time correlators """
    Correlators_equalTime = Matrix{Rational{Int64}}(undef, length(Correlators),max_order+1)
    for j in eachindex(Correlators)
        Correlators_equalTime[j,:] = [sum(Correlators[j][n+1,:] .* [1//1,1//12,1//720,1//30240,1//1209600,1//47900160,691//1307674368000,1//74724249600,3617//10670622842880000,43867//5109094217170944000]) for n in 0:max_order]
    end
    return Correlators_equalTime
end

########## Brillouin Zone functions
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
    """computes the fourier transform along a 1D path through k-space
    given the path computet by create_brillouin_zone_path """
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
 


