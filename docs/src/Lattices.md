## Implementing New Lattice Geometries

This library supports several common lattice geometries out-of-the-box (chain, square, triangular, honeycomb, kagome, pyrochlore, etc.), but you can also implement your own custom lattice geometry. This functionality is built on top of the `Lattice` module from SpinMC.jl.

### How to Add a New Lattice Geometry

To implement a new lattice geometry, you need to add a new case to the `get_finite_Lattice` function in `LatticeGraphs.jl`. Follow these steps:

1. **Define lattice vectors**: These are the primitive vectors of your unit cell
2. **Create a unit cell**: Use the lattice vectors to create a `UnitCell` object
3. **Add basis sites**: Add one or more sites to the unit cell using their fractional coordinates
4. **Define interactions**: Connect sites together to form the lattice structure
5. **Set lattice dimensions**: Define the size of the lattice

### Example

Here's a simple example showing how the square lattice is implemented:

```julia
elseif geometry == "square" ### Square lattice
    a1 = (1, 0)              # First lattice vector
    a2 = (0, 1)              # Second lattice vector
    uc = UnitCell(a1,a2)     # Create the unit cell

    b0 = addBasisSite!(uc, (0.0, 0.0))  # Add a basis site at the origin

    # Add nearest-neighbor interactions
    addInteraction!(uc, b0, b0, (1, 0))  # Connect to neighbor in a1 direction
    addInteraction!(uc, b0, b0, (0, 1))  # Connect to neighbor in a2 direction

    l = (L, L)  # L×L lattice
```

### Implementing a Custom Lattice

To add your own lattice (for example, a "custom_lattice"), add a new `elseif` condition:

```julia
elseif geometry == "custom_lattice"
    # 1. Define lattice vectors
    a1 = (x1, y1)  # Replace with your vector components
    a2 = (x2, y2)  # Replace with your vector components
    uc = UnitCell(a1, a2)
    
    # 2. Add basis sites
    b0 = addBasisSite!(uc, (pos_x0, pos_y0))
    b1 = addBasisSite!(uc, (pos_x1, pos_y1))
    # Add more sites as needed
    
    # 3. Add interactions between sites
    # Format: addInteraction!(uc, source_site, target_site, (unit_cell_offset))
    addInteraction!(uc, b0, b1, (0, 0))  # Connect sites in same unit cell
    addInteraction!(uc, b0, b0, (1, 0))  # Connect to next unit cell in a1 direction
    
    # 4. Set lattice dimensions
    l = (L, L)  # Usually L×L for 2D lattices
```

### Notes

- For 3D lattices, add a third lattice vector `a3` and create the unit cell with `UnitCell(a1, a2, a3)`
- The `addInteraction!` function creates bidirectional connections between sites
- The unit_cell_offset tuples `(n1, n2)` or `(n1, n2, n3)` specify how many unit cells to move in each lattice vector direction

After implementing your lattice, you can create it using:

```julia
mylattice, mygraph = get_finite_Lattice(L, "custom_lattice", PBC=true)
```

## Adding Symmetry Groups for New Lattices

The library includes functionality to reduce computational complexity by exploiting lattice symmetries. This tutorial explains how to define symmetry groups for custom lattices in the `getSymmetryGroup` function.

### Understanding Symmetry Elements

At the core of lattice symmetries are symmetry elements represented by the `sym_element` struct, which consists of two components:

```julia
mutable struct sym_element
    gMat    ::Matrix{Float64}  # Transformation matrix (rotation, reflection, etc.)
    gVec    ::Vector{Float64}  # Translation vector
end
```

When a symmetry element `g` acts on a position vector `r`, it applies:
- First, the rotation/reflection: `gMat * r`
- Then, the translation: `gMat * r + gVec`

For example:
- A pure rotation around the origin has `gVec = [0, 0, ...]`
- A pure translation has `gMat = I` (identity matrix)
- A non-symmorphic symmetry combines both non-trivial `gMat` and `gVec`

### How to Add a New Symmetry Group

To implement a symmetry group for a new lattice, add a new case to the `getSymmetryGroup` function:

1. **Define the lattice vectors** to create the translation group
2. **Create basic symmetry elements** (rotations, reflections, etc.)
3. **Build a basis** of symmetry elements
4. **Generate the full symmetry group**

### Example: Square Lattice Symmetry Group

Here's how the square lattice symmetry group is implemented:

```julia
elseif geometry == "square"
    # 1. Define lattice vectors
    a1 = (1, 0)
    a2 = (0, 1)

    # 2. Create basic symmetry elements
    C_4 = sym_element([0 1; -1 0], [0,0])  # 90° rotation
    Px = sym_element([-1 0; 0 1], [0,0])   # Mirror reflection across y-axis

    # 3. Build a basis with these elements
    basis = sym_group([neutral_elem(C_4), C_4, Px])
    
    # 4. Define the translation group
    translation_Group = translation_group([a1,a2])
    
    # 5. Generate the full symmetry group
    symmetry_Group = generate_symmetry_group(basis, translation_Group)
```

### Adding Your Custom Symmetry Group

To add a symmetry group for your custom lattice:

```julia
elseif geometry == "custom_lattice"
    # 1. Define lattice vectors
    a1 = (...)
    a2 = (...)
    # Add a3 for 3D lattices
    
    # 2. Create symmetry elements
    # Point group operations (rotations, reflections)
    R = sym_element([...], [...])   # Rotation matrix and translation vector
    M = sym_element([...], [...])   # Mirror/reflection
    
    # For non-symmorphic symmetries:
    NS = sym_element([...], [...])  # Matrix and fractional translation
    
    # 3. Create the basis
    basis = sym_group([neutral_elem(R), R, M, NS])
    
    # 4. Define translation group
    translation_Group = translation_group([a1, a2])
    
    # 5. Generate full symmetry group
    symmetry_Group = generate_symmetry_group(basis, translation_Group)
```

## Functions

### Symmetry Elements and Groups

```@docs
DynHTE.sym_element
DynHTE.sym_group
DynHTE.translation_group
DynHTE.neutral_elem
DynHTE.is_element
DynHTE.floor_tol
DynHTE.mod!
DynHTE.add_element!
DynHTE.inverse
DynHTE.commutator
DynHTE.find_order
DynHTE.generate_closed_basis
DynHTE.generate_symmetry_group
DynHTE.shiftRotation
```

### Bonds and Lattice Operations

```@docs
DynHTE.bond
DynHTE.flip_bond!
DynHTE.is_approximately_integer_vector
DynHTE.bond_matrix
```

### Predefined Symmetry Groups and Helpers

```@docs
DynHTE.sym_reduced_lattice
DynHTE.getSymmetryGroup
```


### Lattice Operations

```@docs
DynHTE.Lattice
DynHTE.UnitCell
DynHTE.get_finite_Lattice
DynHTE.latticeToGraph
DynHTE.find_graph_center
DynHTE.find_site_basis_label
```
