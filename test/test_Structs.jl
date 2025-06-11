using SimpleWeightedGraphs
@testset "Structs Tests" begin

    @testset "UnitCell" begin
        # Test 1D constructor
        uc1 = UnitCell((1.0,))
        @test length(uc1.primitive) == 1
        @test length(uc1.basis) == 0
        
        # Test 2D constructor
        uc2 = UnitCell((1.0,0.0), (0.0,1.0))
        @test length(uc2.primitive) == 2
        
        # Test 3D constructor
        uc3 = UnitCell((1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0))
        @test length(uc3.primitive) == 3
        
        # Test type parameters
        @test typeof(uc1) <: UnitCell{1}
        @test typeof(uc2) <: UnitCell{2}
        @test typeof(uc3) <: UnitCell{3}
    end

    @testset "Graph Structures" begin
        # Test Graph construction
        g = SimpleWeightedGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)
        @test_nowarn DynHTE.Graph(g)
        
        # Test GraphG
        @test_nowarn GraphG(g, [1,2])
        
        # Test gG_properties
        prop = gG_properties([1,1,2,0])
        @test prop.order == 1
        @test prop.symmetry_factor == 2
        @test prop.is_symmetric == false
        
        # Test unique_Graph
        graphg = GraphG(g, [1,2])
        @test_nowarn unique_Graph([graphg, [[1,1,2,0]], 1])
        
        # Test unique_Graphs
        @test_nowarn unique_Graphs(2, [])
        
        # Test Dyn_HTE_Graphs
        @test_nowarn Dyn_HTE_Graphs(1//2, unique_Graphs(2,[]), Vector{Vector{Vector{Rational{Int128}}}}())
    end

end


