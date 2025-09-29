# Convenience Functions

This page documents the utility functions provided in the `ConvenienceFunctions.jl` module, which offers various tools for  high-temperature series expansion, correlation function analysis, and spectral function calculations and post-processing.


## Loading Functions

```@docs
getLattice
load_dyn_hte_graphs
Dyn_HTE_Lattice
Dyn_HTE_Graphs
```

## Correlation Functions

```@docs
get_c_iipDyn_mat
get_c_iipEqualTime_mat
```


## Brillouin Zone Path Generation

```@docs
create_brillouin_zone_path
```

## Fourier Transforms

```@docs
get_c_k
inverse_fourier_transform
```

## Sublattice-Resolved Fourier Transforms

```@docs
get_c_k_subl
inverse_fourier_transform_subl
```


## Matsubara Frequency and Polynomial Functions

```@docs
get_TGiip_Matsubara_xpoly
flipEvenIndexEntries
```

## Resummation Tools

```@docs
get_pade
get_p_u
get_LinearTrafoToCoeffs_u
extrapolate_series
```


## Moments and Continued Fractions

```@docs
get_moments_from_c_kDyn
fromMomentsToδ
contFrac
extrapolate_δvec
get_extrapolation_params
```

## Dynamic Structure Factors

```@docs
get_JSkw_mat
ContFracTerminator
JSwithTerminator
```


## Utility Functions

```@docs
create_spin_string
```
