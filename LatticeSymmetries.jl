######## In this file we define symmetry groups to reduce the number of Susceptibilities we have to calculate

using LinearAlgebra

"""
    sym_element

Represents a symmetry element `g`, defined by a linear transformation (matrix `gMat`)
and a translation vector `gVec`. The action of `g` on a vector `r` is:

    g * r = gMat * r + gVec
"""
mutable struct sym_element
    gMat    ::Matrix{Float64}
    gVec    ::Vector{Float64}
end

"""
    sym_group

Represents a symmetry group consisting of a collection of symmetry elements.
"""
mutable struct sym_group
    n_elements   ::Int                 # Number of symmetry elements in the group
    elements     ::Vector{sym_element} # Vector of symmetry elements

    function sym_group(data::Vector{sym_element})
        return new(length(data), data)
    end
end

"""
    translation_group

Represents a group of translational symmetries defined by a basis.
Each column of `basis` is a translation vector (i.e., [a1 a2 ... an]).
"""
mutable struct translation_group
    basis   ::Matrix{Float64}  # Columns are translation basis vectors

    function translation_group(vec_of_tuples::Vector{<:Tuple})
       return new(hcat(collect.(vec_of_tuples)...))
    end
end



"""
    is_element(g, G) -> Bool

Check if symmetry element `g` is part of symmetry group `G`.
"""
function is_element(g::sym_element, G::sym_group)::Bool
    for g2 in G.elements
        if g == g2
            return true
        end
    end
    return false
end

"""version of modf(x) but with a tolarance to counter floating point problems"""
function integer_part_tol(x::Float64; tol=1e-8)
    int_part = floor(Int, x)
    if x - int_part > 1 - tol
        return int_part + 1
    elseif x - int_part < tol - 1
        return int_part - 1
    else
        return int_part
    end
end


"""
    mod!(g, T)

Modifies `g.gVec` so it's within the unit cell defined by translation group `T`.
This ensures the vector is reduced modulo the translation lattice.
"""
function mod!(g::sym_element, T::translation_group)
    nvec = integer_part_tol.(inv(T.basis) * g.gVec) # integer part
    g.gVec .= g.gVec - T.basis * nvec   
    return nothing
end

# Overloading equality for sym_element
import Base.==
function ==(g1::sym_element, g2::sym_element)::Bool
    return isapprox(g1.gVec, g2.gVec; atol = 1e-12)&& isapprox(g1.gMat,g2.gMat; atol = 1e-12)
end

# Overloading equality for sym_group
function ==(G1::sym_group, G2::sym_group)::Bool
    if G1.n_elements == G2.n_elements
        if all(x -> is_element(x, G2), G1.elements)
            return true
        end
    end
    return false
end



# Overloading exponentiation: g^n means applying g n times
import Base.^
function ^(g::sym_element, n::Int)::sym_element
    gout = g
    for i = 1:(n-1)
        gout = g ∘ gout  # Composition
    end
    if n == 0
        gout = neutral_elem(g)
    end
    return gout
end

"""
    neutral_elem(g) -> sym_element

Returns the identity element of the same dimension as `g`.
"""
function neutral_elem(g::sym_element)::sym_element
    d = length(g.gVec)
    sym_element(Matrix(I, d, d), zeros(Float64, d))
end

"""
    neutral_elem(G) -> sym_element

Returns the identity element for the group `G`, based on its first element.
"""
function neutral_elem(G::sym_group)::sym_element
    return neutral_elem(G.elements[1])
end

"""
    g2 ∘ g1 -> sym_element

Group operation: composition of symmetry elements (matrix multiplication and vector addition).
"""
function ∘(g2::sym_element, g1::sym_element)::sym_element
    gMat = g2.gMat * g1.gMat
    gVec = g2.gMat * g1.gVec + g2.gVec
    sym_element(gMat, gVec)
end

"""
    add_element!(g, G)

Adds symmetry element `g` to group `G` and updates its element count.
"""
function add_element!(g::sym_element, G::sym_group)
    push!(G.elements, g)
    G.n_elements += 1
end

"""
    inverse(g) -> sym_element

Returns the inverse of symmetry element `g`.
"""
function inverse(g1::sym_element)
    gMat = inv(g1.gMat)
    gVec = -inv(g1.gMat) * g1.gVec
    sym_element(gMat, gVec)
end

"""
    commutator(g1, g2) -> sym_element

Computes the commutator of two symmetry elements:
    g3 = g2⁻¹ ∘ g1⁻¹ ∘ g2 ∘ g1
"""
function commutator(g1::sym_element, g2::sym_element)
    g1_inv = inverse(g1)
    g2_inv = inverse(g2)
    g3 = g2_inv ∘ g1_inv ∘ g2 ∘ g1
    return g3
end

"""
    find_order(g) -> Int

Finds the smallest positive integer `n` such that g^n = identity.
Throws an error if order > 20.
"""
function find_order(g::sym_element,T::translation_group)::Int
    e = neutral_elem(g)
    order = 0
    gn = e

    while true
        gn = g ∘ gn
        order += 1
        mod!(gn,T)
        if  gn == neutral_elem(g)
            return order
        end
        if order > 20
            throw(error("order of g > 20"))
        end
    end
end

"""
    generate_closed_basis(basis, T) -> sym_group

Generates an overcomplete basis of symmetry elements that is closed under commutation.
Ensures all commutators between elements are included in the group.
"""
function generate_closed_basis(basis::sym_group, T::translation_group)::sym_group
    Gnew = sym_group(copy(basis.elements))

    while true
        G = Gnew
        for (n, g) in enumerate(G.elements)
            for g2 in G.elements[n:end]
                h = commutator(g, g2)
                mod!(h, T)

                if !(is_element(h, G))
                    add_element!(h, Gnew)
                end
            end
        end

        if Gnew == G
            return G
        end

        if Gnew.n_elements > 100
            throw(error("group too large (>100 elements)"))
        end
    end
end

"""
    generate_symmetry_group(basis, T) -> sym_group

Generates the full symmetry group from a given basis of symmetry elements.
This includes elements generated by composition and powers of basis elements.
"""
function generate_symmetry_group(basis::sym_group, T::translation_group)::sym_group
    # First compute a commutation-closed basis
    for element in basis.elements
        mod!(element, T)
    end

    closed_basis = generate_closed_basis(basis, T)
    Gnew = sym_group(copy(closed_basis.elements))

    for g in closed_basis.elements
        G = Gnew
        order = find_order(g,T)
        for n in 0:order
            for g2 in G.elements
                gnew = g^n ∘ g2
                mod!(gnew, T)
                if !(is_element(gnew, G))
                    add_element!(gnew, Gnew)
                end
            end
        end
    end

    return Gnew
end
############Symmetry reduced Lattices

"""
    bond

Represents a bond between two points `r1` and `r2` in the lattice.
"""
mutable struct bond
    r1     ::Vector{Float64}
    r2     ::Vector{Float64}
end

"""
    ==(bond1, bond2) -> Bool

Checks whether two bonds are equal. Two bonds are considered equal if their endpoints
are the same (in either order).
"""
function ==(bond1::bond, bond2::bond)::Bool
    if isapprox(bond1.r1, bond2.r1; atol = 1e-12) && isapprox(bond1.r2, bond2.r2; atol = 1e-12)
        return true
    end
    if isapprox(bond1.r1 , bond2.r2; atol = 1e-12) && isapprox(bond1.r2, bond2.r1; atol = 1e-12)
        return true
    end
    return false
end


""" flip the bond from [r1,r2] to [r2,r1] """
function flip_bond!(b::bond)
    r1 =  b.r1
    r2 =  b.r2
    b.r1 = r2
    b.r2 = r1
  end



"""
    g ∘ b -> bond

Applies a symmetry operation `g` to a bond `b`, transforming both endpoints.
"""
function ∘(g::sym_element, b::bond)::bond
    r1new = g.gMat * b.r1 + g.gVec
    r2new = g.gMat * b.r2 + g.gVec
    bond(r1new, r2new)
end

"""
    mod!(b, T)

Translates a bond `b` so that the first endpoint `r1` lies in the first unit cell
of the translation group `T`. 
"""
function mod!(b::bond, T::translation_group,unitcell::UnitCell)
    nvec = b.r1
    for bas in unitcell.basis
        nvec = inv(T.basis) * (b.r1-collect(bas)) #get integer part  
        if is_approximately_integer_vector(nvec)
            break
        end
    end
    b.r1 .= b.r1 - T.basis * nvec     # Shift r1 into first unit cell
    b.r2 .= b.r2 - T.basis * nvec     # Shift r2 by the same amount to preserve relative position
end


function is_approximately_integer_vector(vec::AbstractVector{<:Float64}; tol=1e-8)
    return all(x -> isapprox(x, round(x); atol=tol), vec)
end

"""
    bond_matrix(lattice, center_sites) -> Matrix{bond}

Constructs a matrix of all bonds between a set of central sites and all sites in the lattice.

Arguments:
- `lattice`: an object representing the full lattice, with a `sitePositions` field.
- `center_sites`: indices of central reference sites.

Returns:
- A matrix where each row corresponds to a central site, and each column to a site in the lattice.
"""
function bond_matrix(lattice::Lattice, center_sites)
    mat = Matrix{bond}(undef, length(center_sites), lattice.length)

    for b in eachindex(center_sites), j in 1:lattice.length
        r1 = collect(lattice.sitePositions[center_sites[b]])
        r2 = collect(lattice.sitePositions[j])
        mat[b, j] = bond(r1, r2)
    end

    return mat
end


"""
    sym_reduced_lattice(lattice, center_sites, G, T)
        -> (reduction_dict, bond_vec_red, position_dict)

Reduces the full list of bonds in the lattice using the provided symmetry group `G`.
Identifies which bonds are equivalent under symmetry operations.

Arguments:
- `lattice`: the full lattice object with `sitePositions`.
- `center_sites`: indices of the reference sites to form initial bonds.
- `G`: symmetry group to apply.
- `T`: translation group to normalize coordinates.

Returns:
- `reduction_dict`: a dictionary mapping each full bond to its reduced (representative) bond index.
- `bond_vec_red`: list of symmetry-inequivalent bonds (the reduced basis).
- `position_dict`: mapping from reduced bond index to lattice site indices (origin, target).
"""
function sym_reduced_lattice(lattice::Lattice, center_sites, G::sym_group, T::translation_group)
    bmat = bond_matrix(lattice, center_sites)             # Matrix of all bonds
    bond_vec_red = copy(bmat[:])                          # Flattened list of all bonds
    reduction_dict = Dict()                               # Map from full bond index to reduced bond index
    position_dict = Dict()                                # Map from reduced bond index to actual site indices

    #remove duplicates inside the unit cell

    for (it, b) in enumerate(bond_vec_red)
        # Find original (row, col) location of bond in matrix
        orig_ind = findfirst(x -> x == b, bmat)
        push!(position_dict, it => (center_sites[orig_ind[1]], orig_ind[2]))
        

        for g in G.elements
            #find new element
            bnew = g ∘ b
            mod!(bnew, T,lattice.unitcell)  # Normalize to unit cell
            
            # check its index
            mat_ind = findall(x -> x == bnew, bmat)

           
            # Record symmetry equivalence
            for ind in  mat_ind
            push!(reduction_dict, ind => it)
            end

            if  b≠ bnew
                #remove the corresponding bonds from bond_vec_red, that have a higher index than the original bond
                todelete = findall(x -> x == bnew, bond_vec_red)
                filter!( x -> x > it,todelete)
                deleteat!(bond_vec_red,todelete)
            end

            #do the same for the fliped bond
            flip_bond!(bnew)

            mod!(bnew, T,lattice.unitcell)

            mat_ind = findall(x -> x == bnew, bmat)

            for ind in  mat_ind
            push!(reduction_dict, ind => it)
            end

            if  b≠ bnew
                todelete = findall(x -> x == bnew, bond_vec_red)
                filter!( x -> x > it,todelete)
                deleteat!(bond_vec_red,todelete)
            end

        end
        
    end

    return reduction_dict, bond_vec_red, position_dict
end

function getSymmetryGroup(geometry::String)
    """ creates symmetry group for corresponding Lattice
    Currently implmented:
        - square
        - triang
        - honeycomb
        - pyrochlore
        - kagome
    """
    if geometry == "chain" ### chain lattice
        a1 = (1, 0)
        a2 = (0, 1)

        Px = sym_element([-1 0; 0 1],[0,0])
        basis = sym_group([neutral_elem(Px),Px])
        translation_Group = translation_group([a1,a2])
        symmetry_Group = generate_symmetry_group(basis,translation_Group)

    elseif geometry == "square" ### Square lattice
        a1 = (1, 0)
        a2 = (0, 1)

        C_4 = sym_element([0 1; -1 0],[0,0])
        Px = sym_element([-1 0; 0 1],[0,0])

        basis = sym_group([neutral_elem(C_4),C_4,Px])
        translation_Group = translation_group([a1,a2])
        symmetry_Group = generate_symmetry_group(basis,translation_Group)

    elseif geometry == "simple_cubic" ### Square lattice
        a1 = (1, 0, 0)
        a2 = (0, 1 , 0)
        a3 = (0, 0 , 1)

        C_41 = sym_element([0 1 0; -1 0 0; 0 0 1],[0,0,0])
        C_42 = sym_element([0 0 1; 0 1 0; -1 0 0],[0,0,0])
        Px = sym_element([-1 0 0 ; 0 1 0; 0 0 1],[0,0,0])

        basis = sym_group([neutral_elem(C_41),C_41,C_42,Px])
        translation_Group = translation_group([a1,a2,a3])
        symmetry_Group = generate_symmetry_group(basis,translation_Group)

    elseif geometry == "triang" ### Square lattice
        a1 = (1/2, sqrt(3)/2)
        a2 = (1/2, -sqrt(3)/2)

        C_6 = sym_element([1/2 -sqrt(3)/2; sqrt(3)/2 1/2],[0,0])
        Px = sym_element([-1 0; 0 1],[0,0])
        basis = sym_group([neutral_elem(C_6),C_6,Px])
        translation_Group = translation_group([a1,a2])
        symmetry_Group = generate_symmetry_group(basis,translation_Group)

    elseif geometry == "kagome" ### Square lattice
        a1 = (1, sqrt(3))
        a2 = (1, -sqrt(3))

        C_6 = shiftRotation([1/2 -sqrt(3)/2; sqrt(3)/2 1/2],[-1,0]) #C6 rotation around center of hexagon in Kagome lattice. 
        Px = sym_element([-1 0; 0 1],[0,0])
        basis = sym_group([neutral_elem(C_6),C_6,Px])
        translation_Group = translation_group([a1,a2])
        symmetry_Group = generate_symmetry_group(basis,translation_Group)

    elseif geometry == "pyrochlore" ### Square lattice
        a1 = (0,1/2,1/2)
        a2 = (1/2,0,1/2)
        a3 = (1/2,1/2,0)
        translation_Group = translation_group([a1,a2,a3])

        C_6I = sym_element([0 -1 0 ; 0 0 -1 ; -1 0 0],[0,0,0]) #C6 rotoreflection. 
        S = sym_element([0 0 1 ; 0 -1 0 ; 1 0 0],[1/4,0,1/4])#Screw symmetry along (1,0,1) 
        basis = sym_group([neutral_elem(C_6I),C_6I,S])
        symmetry_Group = generate_symmetry_group(basis,translation_Group)

    
    else 
        throw(error("Symmetry group for geometry: " * geometry * " not yet implemented"))
    end

    ### create graph for the lattice and return
    return symmetry_Group,translation_Group
end

""" takes a rotation R about a point p and outputs the symmetry element that rotates around the origin and translates by p-Rp """
function shiftRotation(R::Matrix{<:Number},p::Vector{<:Number})::sym_element
   return  sym_element(R,p-R*p)
end


#= 

sym_G,transl_G = getSymmetryGroup("pyrochlore")

@time reduction_dict,bond_vec_red,position_dict = sym_reduced_lattice(hte_lattice.lattice,hte_lattice.basis_positions,sym_G,transl_G)

reduction_dict
bond_vec_red
position_dict
 =#

#= 
sym_G,transl_G = getSymmetryGroup("triang")
sym_G,transl_G = getSymmetryGroup("square")

reduction_dict,bond_vec_red,position_dict = sym_reduced_lattice(hte_lattice.lattice,center_sites,sym_G,transl_G)


reduction_dict



function compute_lattice_correlations(LatGraph,lattice,reduction_dict,bond_vec_red,position_dict,max_order,gG_vec_unique,C_Dict_vec)::Array{Matrix{Rational{Int128}}}
    """compute all correlations from the center_sites to all other sites of the lattice"""
    Correlators = Array{Matrix{Rational{Int128}}}(undef, lattice.length,length(lattice.unitcell.basis));
    reduced_Correlators = Vector{Matrix{Rational{Int128}}}(undef, length(bond_vec_red));

    Threads.@threads for i in reverse(eachindex(bond_vec_red))
        reduced_Correlators[i] = mapreduce(permutedims, vcat, Calculate_Correlator_fast(LatGraph,position_dict[i][1],position_dict[i][2],max_order,gG_vec_unique,C_Dict_vec))
    end
    
    for jp = 1:lattice.length
        for b = 1:length(lattice.unitcell.basis)
        Correlators[jp,b] = reduced_Correlators[reduction_dict[CartesianIndex(b,jp)]]
        end
    end
    return Correlators
end


@time Correlators1 = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);
@time Correlators2 = compute_lattice_correlations(LatGraph,lattice,reduction_dict,bond_vec_red,position_dict,max_order,gG_vec_unique,C_Dict_vec);
Correlators1 == Correlators2 =#