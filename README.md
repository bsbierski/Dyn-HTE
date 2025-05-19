# Dyn-HTE (post-processing)
High temperature expansion of dynamic Matsubara spin correlators for general Heisenberg models (only post-processing)

* **Dyn-HTE:** The main code which shows a usage example of Dyn-HTE both for real- and momentum-space Matsubara correlators

* **LatticeGraphs** Generates lattices or lattice balls L by keyword geometry=chain,square,triangle,... (can be expanded by user). The site indices of the center unit-cell of the lattice ball are returned. This code builds on Lattice.jl which is taken from SpinMC.jl

* **Embedding:**  For given L, calculation of embedding factors for graphG and calculation of expansion coefficients of $TG_{ii\prime}(i\nu_m)$ (expansion in powers of -x)

* **ConvenienceFunctions:**  Definition of various functions to help evaluate and plot the results of the Dyn-HSE

* **Dyn-HTE_Tutorial_TriangularLattice.pdf** Tutorial on how to use the Dyn-HTE working through the example of the S=1/2 nearest-neighbor Heisenberg AFM on the triangular lattice.