
```@meta
CurrentModule = DynHTE
```


# DynHTE

Documentation for [DynHTE](https://github.com/bsbierski/Dyn-HTE.jl).


# Dynamic high-temperature expansion for quantum spins (Dyn-HTE) 

by Ruben Burkard, Benedikt Schneider, Björn Sbierski 

This software allows to compute dynamic spin correlations of spin-S Heisenberg models in thermal equilibrium via a high-temperature expansion (HTE) of the Matsubara spin-spin correlator. It is currently applicable for ``S \leq 1`` models with a single coupling constant ``J``, but arbitrary lattice geometry.


```math
H = J \sum_{(ii^\prime)}\left(S_{i}^{+}S_{i^\prime}^{-} + S_{i}^{-}S_{i^\prime}^{+} + S_{i}^{z}S_{i^\prime}^{z}\right)
```

Dyn-HTE obtains the HTE of the Matsubara spin-spin correlator

```math
G_{ii^\prime}(i\nu_{m}) = T \int_{0}^{\beta} \mathrm{d}\tau\,\mathrm{d}\tau^{\prime}\,e^{i\nu_{m}(\tau - \tau^{\prime})}
\left\langle \mathcal{T} S_{i}^{z}(\tau) S_{i^{\prime}}^{z}(\tau^{\prime}) \right\rangle
```

and allows post-processing of this information to compute the dynamic spin structure factor (DSF):

```math
S(\mathbf{k}, \omega) = \int_{-\infty}^{+\infty} \frac{\mathrm{d}t}{2\pi N} \sum_{i,i^\prime} 
\exp(i\omega t - i\mathbf{k} \cdot (r_i - r_{i^\prime})) \left\langle S_{i}^{z}(t) S_{i^\prime}^{z} \right\rangle
```



## Installation

The DynHTE.jl package can be installed by invoking the following command in the Julia REPL:
```julia
]  add https://github.com/bsbierski/Dyn-HTE/tree/package
```
Due to the size of the precomputed data, this can take a minute. 


## Publication/Citation
The theory background for Dyn-HTE and various applications are provided in the following two publications:

[1] Ruben Burkard, Benedikt Schneider, Björn Sbierski, *High-temperature series expansion of the dynamic Matsubara spin correlator*, arxiv 2505.23699 (2025)

[2] Ruben Burkard, Benedikt Schneider, Björn Sbierski, *Dynamic correlations of frustrated quantum spins from high-temperature expansion*, arxiv 2505.14571 (2025)

If Dyn-HTE benefits your research, please acknowledge it by citing these references.

# Index
```@index
```