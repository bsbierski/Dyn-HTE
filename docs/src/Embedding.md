# Graph Embedding and Generation

This page documents the functions used for graph generation, manipulation, and embedding calculations in the high-temperature series expansion.

## Embedding of Correlators

```@docs
DynHTE.Calculate_Correlator_fast
DynHTE.e_fast
```

## Graph Structs

```@docs
DynHTE.unique_Graph
DynHTE.unique_Graphs
DynHTE.GraphG
DynHTE.Graph
DynHTE.gG_properties
```

## Graph Isomorphism and Structure
```@docs
DynHTE.is_simple_isomorphic
DynHTE.isIsomorph
DynHTE.is_symmetric
DynHTE.give_unique_gG_vec
```

## Graph Properties
```@docs
DynHTE.degeneracy
DynHTE.totalEdges
DynHTE.noLeavesExceptAt
DynHTE.hasGeneralizedLeaves
DynHTE.symmetryFactor
```

## Graph Transformation and Manipulation

```@docs
DynHTE.toSimpleGraph
```

## Graph Generation and Collection

```@docs
DynHTE.getVacGraphs
DynHTE.getGraphsG
```

## Graph Visualization


```@docs
DynHTE.gplot
```
