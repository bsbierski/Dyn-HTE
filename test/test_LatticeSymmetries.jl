

@testset "LatticeSymmetries.jl Tests" begin
    
    @testset "Symmetry Element Operations" begin
        # Create basic symmetry elements
        id = sym_element(Matrix{Float64}(I, 2, 2), [0.0, 0.0])
        reflection = sym_element([-1.0 0.0; 0.0 1.0], [0.0, 0.0])
        rotation = sym_element([0.0 -1.0; 1.0 0.0], [0.0, 0.0])  # 90-degree rotation
        translation = sym_element(Matrix{Float64}(I, 2, 2), [1.0, 0.0])
        
        # Test equality
        @test id == sym_element(Matrix{Float64}(I, 2, 2), [0.0, 0.0])
        @test reflection != rotation
        
        # Test composition
        @test reflection ∘ reflection == id
        @test rotation ∘ rotation ∘ rotation ∘ rotation == id
        @test translation ∘ translation == sym_element(Matrix{Float64}(I, 2, 2), [2.0, 0.0])
        
        # Test inverse
        @test inverse(id) == id
        @test inverse(reflection) == reflection
        @test inverse(rotation) == sym_element([0.0 1.0; -1.0 0.0], [0.0, 0.0])
        @test inverse(translation) == sym_element(Matrix{Float64}(I, 2, 2), [-1.0, 0.0])
        
        # Test commutator
        @test commutator(id, reflection) == id
        @test commutator(rotation, reflection) != id  # These don't commute
        
        # Test power operation
        @test reflection^2 == id
        @test rotation^4 == id
        @test rotation^2 == sym_element([-1.0 0.0; 0.0 -1.0], [0.0, 0.0])
        
        # Test neutral element
        @test neutral_elem(rotation) == id
    end
    
    @testset "Symmetry Group Operations" begin
        id = sym_element(Matrix{Float64}(I, 2, 2), [0.0, 0.0])
        reflection = sym_element([-1.0 0.0; 0.0 1.0], [0.0, 0.0])
        
        # Create a simple Z2 group
        z2_group = sym_group([id, reflection])
        
        # Test element inclusion
        @test is_element(id, z2_group)
        @test is_element(reflection, z2_group)
        @test !is_element(sym_element([0.0 -1.0; 1.0 0.0], [0.0, 0.0]), z2_group)
        
        # Test group equality
        @test z2_group == sym_group([id, reflection])
        @test z2_group != sym_group([id])
        
        # Test neutral element
        @test neutral_elem(z2_group) == id
        
        # Test adding elements
        g = sym_group([id])
        add_element!(reflection, g)
        @test g == z2_group
        @test g.n_elements == 2
    end
    
    @testset "Translation Group" begin
        # Create a 2D square lattice translation group
        square_trans = translation_group([(1.0, 0.0), (0.0, 1.0)])
        
        # Test that basis is correctly set up
        @test size(square_trans.basis) == (2, 2)
        @test square_trans.basis[:,1] ≈ [1.0, 0.0]
        @test square_trans.basis[:,2] ≈ [0.0, 1.0]
        
        # Test constructor with different tuple types
        int_trans = translation_group([(1, 0), (0, 1)])
        @test int_trans.basis ≈ square_trans.basis
        
        # Test mod! operation for symmetry elements
        g = sym_element(Matrix{Float64}(I, 2, 2), [2.3, 1.7])
        mod!(g, square_trans)
        @test g.gVec ≈ [0.3, 0.7] atol=1e-12
        
        # Test floor_tol function
        @test floor_tol(2.999, tol=0.01) == 3
        @test floor_tol(2.001, tol=0.01) == 2
        @test floor_tol(-2.001, tol=0.01) == -2
        @test floor_tol(-2.999, tol=0.01) == -3
    end
    
    @testset "Group Generation" begin
        # Define translation group for a square lattice
        T = translation_group([(1.0, 0.0), (0.0, 1.0)])
        
        # Define C4 rotation symmetry
        C4 = sym_element([0.0 -1.0; 1.0 0.0], [0.0, 0.0])
        
        # Create a basis with just C4
        basis = sym_group([neutral_elem(C4), C4])
        
        # Generate the full C4 group
        C4_group = generate_symmetry_group(basis, T)
        
        # Should contain 4 elements: 1, C4, C4^2, C4^3
        @test C4_group.n_elements == 4
        
        # Test that all powers of C4 are in the group
        @test is_element(C4^2, C4_group)
        @test is_element(C4^3, C4_group)
        
        # Test finding order
        @test find_order(C4, T) == 4
        @test find_order(C4^2, T) == 2
        @test find_order(C4^3, T) == 4
    end
    
    @testset "Bond Operations" begin
        # Create some bonds
        b1 = bond([0.0, 0.0], [1.0, 0.0])
        b2 = bond([0.0, 0.0], [0.0, 1.0])
        
        # Test bond equality
        @test b1 != b2
        @test b1 == bond([0.0, 0.0], [1.0, 0.0])
        @test b1 == bond([1.0, 0.0], [0.0, 0.0])  # Order doesn't matter
        
        # Test symmetry operation on bonds
        reflection = sym_element([-1.0 0.0; 0.0 1.0], [0.0, 0.0])
        reflected_b1 = reflection ∘ b1
        @test reflected_b1 == bond([0.0, 0.0], [-1.0, 0.0])
        
        # Test bond flipping
        flipped_b1 = bond([1.0, 0.0], [0.0, 0.0])
        flip_bond!(flipped_b1)
        @test flipped_b1 == bond([0.0, 0.0], [1.0, 0.0])
        @test b1 == flipped_b1
    end
    
    @testset "Predefined Symmetry Groups" begin
        # Test that we can generate the symmetry groups for different lattices
        sq_group, sq_trans = getSymmetryGroup("square")
        @test sq_group.n_elements == 8  # D4 has 8 elements
        
        # Triangular lattice has D6 symmetry (12 elements)
        tri_group, tri_trans = getSymmetryGroup("triang")
        @test tri_group.n_elements == 12
        
        # Chain should have just reflection symmetry (2 elements)
        chain_group, chain_trans = getSymmetryGroup("chain")
        @test chain_group.n_elements == 2
        
        # Test simple cubic
        sc_group, sc_trans = getSymmetryGroup("simple_cubic")
        @test sc_group.n_elements > 8  # Should have more elements than square
        
        # Test that an invalid geometry throws an error
        @test_throws ErrorException getSymmetryGroup("invalid_geometry")
    end
    
    @testset "shiftRotation function" begin
        # Test the shiftRotation helper function
        R = [0.0 -1.0; 1.0 0.0]  # 90° rotation
        p = [1.0, 0.0]  # Shift point
        
        shifted_rotation = shiftRotation(R, p)
        
        # The resulting symmetry element should rotate around p
        # To verify: applying it to p should return p
        result = shifted_rotation.gMat * p + shifted_rotation.gVec
        @test isapprox(result, p, atol=1e-12)
        
        # Another point should be rotated around p
        q = [2.0, 0.0]
        expected = [1.0, 1.0]  # q rotated 90° around p
        result = shifted_rotation.gMat * q + shifted_rotation.gVec
        @test isapprox(result, expected, atol=1e-12)
    end
end;
