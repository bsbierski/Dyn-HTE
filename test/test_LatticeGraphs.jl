@testset "LatticeGraphs Tests" begin

    @testset "latticeToGraph Basic Functionality" begin
        # Create a simple lattice (chain) for testing
        lattice, graph = get_finite_Lattice(5, "chain", PBC=false)
        
        # Test conversion to graph
        test_graph = latticeToGraph(lattice)
        
        # Check that the graph has the correct number of vertices
        @test nv(test_graph) == lattice.length
        
        # Check that edges correspond to interactions
        for v in 1:nv(test_graph)
            for neighbor in lattice.interactionSites[v]
                if neighbor != v # Skip self-interactions
                    @test has_edge(test_graph, v, neighbor)
                end
            end
        end
    end

    @testset "Tuple Operations" begin
        # Test subtraction for 2D tuples
        a2d = (3.0, 4.0)
        b2d = (1.0, 2.0)
        @test a2d - b2d == (2.0, 2.0)
        
        # Test subtraction for 3D tuples
        a3d = (3.0, 4.0, 5.0)
        b3d = (1.0, 2.0, 3.0)
        @test a3d - b3d == (2.0, 2.0, 2.0)
    end

    @testset "find_graph_center Function" begin
        # Create a small test graph where the center is obvious
        g = SimpleGraph(5)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 3, 5)
        
        # Node 3 should be the center
        @test 3 in find_graph_center(g)
        
        # Create another graph - a path graph where middle node is center
        path = path_graph(5)
        @test 3 in find_graph_center(path)

        # Test complete graphs where all nodes are centers
        for n = 3:5
            g = complete_graph(n)
            @test find_graph_center(g) == collect(1:n)
        end
    end

    @testset "getLattice Function Basics" begin
        # Test for various geometries
        for geometry in ["chain", "square", "triang"]
            L = 3 # Small L for testing
            dyn_hte_lattice = getLattice(L, geometry)
            
            # Check structure correctness
            @test dyn_hte_lattice.name == geometry
            @test isa(dyn_hte_lattice.lattice, Lattice)
            @test isa(dyn_hte_lattice.graph, SimpleGraph)
            @test !isempty(dyn_hte_lattice.basis_positions)
            
            # Check graph properties
            @test nv(dyn_hte_lattice.graph) == dyn_hte_lattice.lattice.length
            @test is_connected(dyn_hte_lattice.graph)
        end
    end

    @testset "Maximum Distance Constraint" begin
        # Test the constraint that from any basis position, distance to any site is at most L
        
        # Test for square lattice
        L = 4
        square_lattice = getLattice(L, "square")
        square_graph = square_lattice.graph
        basis_positions = square_lattice.basis_positions
        
        # Calculate distances from each basis position
        for pos in basis_positions
            distances = dijkstra_shortest_paths(square_graph, pos).dists
            @test all(d -> d <= L, filter(isfinite, distances))
            # Also verify we have at least one site at maximum distance L
            @test L in distances
        end
        
        # Test for triangular lattice
        tri_lattice = getLattice(L, "triang")
        tri_graph = tri_lattice.graph
        tri_basis_positions = tri_lattice.basis_positions
        
        for pos in tri_basis_positions
            distances = dijkstra_shortest_paths(tri_graph, pos).dists
            @test all(d -> d <= L, filter(isfinite, distances))
            @test L in distances
        end
    end

    @testset "find_site_basis_label Function" begin
        L = 3
        # Test for a geometry with multiple basis sites
        honeycomb_lattice = getLattice(L, "honeycomb")
        basis_labels = find_site_basis_label(honeycomb_lattice.lattice)
        
        # Check that all sites have valid basis labels
        @test all(label -> label > 0, basis_labels)
        
        # For a honeycomb lattice, we should have 2 unique basis labels
        @test length(unique(basis_labels)) == 2
    end
end

