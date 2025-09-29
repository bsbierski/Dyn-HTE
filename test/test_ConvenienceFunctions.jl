@testset "Convenience Functions" begin
    @testset "Spin String Creation" begin
        @test create_spin_string(0.5) == "S1half"
        @test create_spin_string(1) == "S1"
        @test create_spin_string(3/2) == "S3half"
        @test create_spin_string(2) == "S2"
        @test_throws ErrorException create_spin_string(1/3)
    end
    
    @testset "Vector Operations" begin
        # Test flipEvenIndexEntries
        @test flipEvenIndexEntries([1, 2, 3, 4]) == [1, -2, 3, -4]
        @test flipEvenIndexEntries([5]) == [5]
        @test flipEvenIndexEntries([1, 2, 3]) == [1, -2, 3]
    end
    
    @testset "Polynomial Transformations" begin
        # Test get_LinearTrafoToCoeffs_u with simple cases
        trafo = get_LinearTrafoToCoeffs_u(2, 1.0)
        @test size(trafo) == (3, 3)
        @test trafo[1,1] ≈ 1.0
        @test trafo[2,2] ≈ 1.0
        
        # Test with different scaling factor
        trafo2 = get_LinearTrafoToCoeffs_u(2, 2.0)
        @test trafo2[2,2] ≈ 0.5 # Should scale by 1/f
        
        # Test maximum order constraint
        @test_throws ErrorException get_LinearTrafoToCoeffs_u(17, 1.0)
    end


    @testset "Continued Fraction Terminator" begin

       F(x, c , σ) = ContFracTerminator(x, σ,c)

        g(x) = begin
             c = BigFloat(4.123)
                σ = BigFloat(2.12512)
            return x * F(x, c,σ) + (c + 1*σ) * F(x, c,σ) * F(x, c + 1*σ,σ) - 1
        end
        xs = range(-5.0, 5.0; length=20)
        ys = round.([real(g(x)) for x in xs]; digits=10)
        @test ys ≈ zeros(length(xs))
    end
    
    @testset "Brillouin Zone Path" begin
        # Test with 2D points
        points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]
        path, indices = create_brillouin_zone_path(points, 10)
        
        @test length(path) ≥ 10
        @test length(indices) == length(points)
        @test path[indices[1]] == points[1]
        @test path[indices[end]] == points[end]
        
        # Test with 3D points
        points3d = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)]
        path3d, indices3d = create_brillouin_zone_path(points3d, 15)
        
        @test length(path3d) ≥ 15
        @test length(indices3d) == length(points3d)
        @test path3d[indices3d[1]] == points3d[1]
        @test path3d[indices3d[end]] == points3d[end]
    end

    
    @testset "Moments and Delta Parameters" begin
        # Create a simple test matrix
        c_kDyn = zeros(13, 7)
        # Set up simple coefficients (identity matrix-like)
        for i in 1:13
            for j in 1:7
                if i == 2*j-1
                    c_kDyn[i,j] = 1.0
                end
            end
        end
        
        # Test moment extraction
        moments = get_moments_from_c_kDyn(c_kDyn)
        @test length(moments) > 0
        @test Polynomials.degree(moments[1]) >= 0
        
        # Test delta conversion with real numbers
        m_vec = [4.0, 2.0, 1.0, 0.5] # Simple test moments
        δ_vec, r_vec = fromMomentsToδ(m_vec)
        @test length(δ_vec) == length(m_vec)
        @test δ_vec[1] == 4.0 # First delta should equal m0
        @test r_vec == [0.0, 1.0, 2.0, 3.0]
        
        # Check the conversion preserves continued fraction representation
        # by testing a simple harmonic oscillator case
    end
    
    
    @testset "Equal Time Correlators" begin
        # Create a simple test matrix
        c_iipDyn_mat = Matrix{Matrix{Rational{Int128}}}(undef, 2, 2)
        
        # Fill with simple test data
        for i in 1:2, j in 1:2
            c_iipDyn_mat[i,j] = ones(Rational{Int128}, 5, 10)
        end
        
        # Test frequency sum
        c_equaltime = get_c_iipEqualTime_mat(c_iipDyn_mat)
        @test size(c_equaltime) == (2, 2)
        @test length(c_equaltime[1,1]) == 5
        
        # The result should be the sum of coefficients multiplied by moment factors
        @test all(x -> x > 0, c_equaltime[1,1])  # All entries should be positive
    end
    
    @testset "Matsubara Correlator" begin
        # Create test data for a simple case
        c_iipDyn_mat = Matrix{Matrix{Rational{Int128}}}(undef, 2, 2)
        
        # Fill with simple test data
        test_matrix = zeros(Rational{Int128}, 12, 10)
        test_matrix[12,1] = 1 // 1
        test_matrix[11,2] = 1 // 1
        
        for i in 1:2, j in 1:2
            c_iipDyn_mat[i,j] = copy(test_matrix)
        end
        
        # Test Matsubara frequency m=0
        poly_m0 = get_TGiip_Matsubara_xpoly(c_iipDyn_mat, 1, 1, 0)
        @test Polynomials.degree(poly_m0) == 11
        @test poly_m0(0.0) == 0.0
        
        # Test non-zero Matsubara frequency
        poly_m1 = get_TGiip_Matsubara_xpoly(c_iipDyn_mat, 1, 1, 1)
        @test Polynomials.degree(poly_m1) == 10
    end
end





       