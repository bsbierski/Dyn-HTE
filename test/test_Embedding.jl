
@testset "Embedding Tests" begin

    @testset "is_simple_isomorphic" begin
        # Create two identical simple graphs
        g1 = SimpleWeightedGraph(3)
        add_edge!(g1, 1, 2)
        add_edge!(g1, 2, 3)
        graphG1 = GraphG(g1, [1, 2])
        graphG2 = GraphG(g1, [1, 2])
        
        @test is_simple_isomorphic(graphG1, graphG2) == true
        
        # Test onsite vs offsite
        g2 = SimpleWeightedGraph(3)
        add_edge!(g2, 1, 2)
        add_edge!(g2, 2, 3)
        graphG3 = GraphG(g2, [1, 1])  # onsite
        graphG4 = GraphG(g2, [1, 2])  # offsite
        
        @test is_simple_isomorphic(graphG3, graphG4) == false
    end

    @testset "e_fast" begin
        n_max = 4
        spin_length = 1
        Lgraph = complete_graph(4)
        hte_graphs = load_dyn_hte_graphs(spin_length,n_max);

        g1 = hte_graphs.unique_graphs.graphs[1].ref_graph
        g2 = hte_graphs.unique_graphs.graphs[2].ref_graph
        g3 = hte_graphs.unique_graphs.graphs[3].ref_graph
        g4 = hte_graphs.unique_graphs.graphs[4].ref_graph
 
        @test e_fast(Lgraph,1,1,g1) == 1
        @test e_fast(Lgraph,1,2,g2) == 1
        @test e_fast(Lgraph,1,1,g3) == 3
        @test e_fast(Lgraph,1,2,g4) == 2
    end

    @testset "Calculate_Correlator_fast" begin
       ### maximum expansion order and spin-length
        n_max = 6
        spin_length = 1/2
        Lgraph = complete_graph(4)
        hte_graphs = load_dyn_hte_graphs(spin_length,n_max);

        c_iipDyn_mat = get_c_iipDyn_mat(Lgraph,[1],hte_graphs)

        uniform_susceptibility = (c_iipDyn_mat[1,1]+c_iipDyn_mat[2,1]+c_iipDyn_mat[3,1]+c_iipDyn_mat[4,1])

        ### test if uniform susceptibility is pureyl static (m=0 only)
        @test unique(uniform_susceptibility[:,2:end]) == Rational{Int128}[0] 
        
         expected_mat_1 = [
            1//4     0        0         0      0  0  0  0  0  0
            0        0        0         0      0  0  0  0  0  0
            -1//32   3//8     0         0      0  0  0  0  0  0
            -1//128  3//32    0         0      0  0  0  0  0  0
            23//2560 -11//128 -21//16   0      0  0  0  0  0  0
            199//30720 -35//512 -9//16  0      0  0  0  0  0  0
            -1859//860160 187//10240 87//256 81//16 0 0 0 0 0 0
        ]
        
         expected_mat_2 = [
            0          0           0          0       0  0  0  0  0  0
            1//16      0           0          0       0  0  0  0  0  0
            5//192    -1//8        0          0       0  0  0  0  0  0
            -3//256   -1//32       0          0       0  0  0  0  0  0
            -221//15360 11//384    7//16      0       0  0  0  0  0  0
            163//184320 35//1536   3//16      0       0  0  0  0  0  0
            11219//1720320 -187//30720 -29//256 -27//16 0 0 0 0 0 0
        ]

        @test c_iipDyn_mat[1,1] == expected_mat_1
        @test c_iipDyn_mat[2,1] == expected_mat_2
    end


end
