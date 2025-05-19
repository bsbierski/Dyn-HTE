# Dynamic high temperature expansion for quantum spins <br /> (Dyn-HTE) 

software by Ruben Burkard, Benedikt Schneider, Björn Sbierski 

This software allows to compute dynamic spin correlations of spin-S Heisenberg models in thermal equilibrium via a high-temperature expansion of the Matsubara spin-spin correlator. It currently allows for $S \leq 1$ and models with a single coupling constant $J$, but arbitrary lattice geometry.

$$ H=J\sum_{(ii^\prime)}\left(S_{i}^{+}S_{i^\prime}^{-}+S_{i}^{-}S_{i^\prime}^{+}+S_{i}^{z}S_{i^\prime}^{z}\right) $$

Dyn-HTE obtains the HTE of the Matsubara spin-spin correlator

$$ G_{ii^\prime}(i\nu_{m})=T  \int_{0}^{\beta} \mathrm{d}\tau\mathrm{d}\tau^{\prime}\,e^{i\nu_{m}(\tau-\tau^{\prime})} 
\left\langle \mathcal{T}S_{i}^{z}(\tau)S_{i^{\prime}}^{z}(\tau^{\prime})\right\rangle $$

and allows to post-process this information into the dynamic spin structure factor (DSF)

$$ S(\mathbf{k},\omega) = \int_{-\infty}^{+\infty}  \frac{\mathrm{d}t}{2\pi N} \sum_{i,i^\prime}  exp(i\omega t-i\mathbf{k}\cdot(r_i - r_{i^\prime}))  \left\langle S_{i}^{z}(t)S_{i^\prime}^{z}\right\rangle $$


## Publication/Citation
The theory background for Dyn-HTE and various applications are provided in the following two publications.

[1] Ruben Burkard, Benedikt Schneider, Björn Sbierski, *Dyn-HTE: High-temperature expansion of the dynamic Matsubara spin correlator*, arxiv 2025.YYYYY (2025)
[2] Ruben Burkard, Benedikt Schneider, Björn Sbierski, *Dynamic correlations of frustrated quantum spins from high-temperature expansion*, arxiv 2025.XXXXX (2025)

If Dyn-HTE benefits your research, please acknowledge it by citing these references.

## Tutorial: Spin-1/2 Heisenberg AFM on triangular lattice
This tutorial explains the use of the Dyn-HTE software provided in this repository using the example of the nearest-neighbor S=1/2 Heisenberg AFM on the triangular lattice. The associated julia script can be found under “CaseStudy/Triangular_Lattice/Application_Triangular_Lattice.jl”. This script contains complete code, here we only highlight the most important functionalities specific to Dyn-HTE and assume the reader is familiar with the julia language and its plotting routines. The physical background and most of the results generated in this tutorial are discussed in the two publications mentioned in [Publication/Citation](#Publication/Citation).

### Preparations: Define lattice and find Dyn-HTE for Matsubara correlator
To start, we need to include the necessary supporting julia scripts and the packages (JLD2, DelimitedFiles) that manage file handling. This is the same for every application of Dyn-HTE.
```bash
using JLD2, DelimitedFiles
include("../../plotConventions.jl")
include("../../LatticeGraphs.jl")
include("../../Embedding.jl")
include("../../ConvenienceFunctions.jl")
```

Next we fix the spin length to $S=1/2$ (also $S=1$ would be currently available) and load all the graph evaluations for this S for the maximum available order n_max=12.
```bash
spin_length = 1/2
n_max = 12
hte_graphs = load_dyn_hte_graphs(spin_length,n_max);
```
The triangular lattice is predefined in the script “LatticeGraphs.jl” and due to the maximum order n_max=12 we only need a piece of it with sites separated from a central site by L=n_max=12 nearest-neighbor bonds or less. Other predefined lattices available via their keyword are “chain”, “square”, “triang”(ular), “honeycomb”, “pyrochlore”, “kagome”. Other translation invariant geometries can be defined by adapting the function “get_finite_Lattice” in “LatticeGraphs.jl”. We also define the three special points in the Brillouin zone (BZ) that will be of interest to us later, the origin $\Gamma$  and $K=(2\pi/3,2\pi/\sqrt{3})$ and $M=(0,2\pi/\sqrt{3})$.

```bash
L = 12
hte_lattice = getLattice(L,"triang");
Γ,K,M = (0,0), (2*π/3,2*π/sqrt(3)), (0,2*π/sqrt(3))
```
Note that it is not necessary to proceed with a Dyn_HTE_Lattice structure, one can also define a SimpleGraph type of the “Graphs.jl” package which is useful for finite or irregular systems (see “SpinCluster_BSb.jl”). 
The lattice and in particular the site numbering can be visualized by the following line (it is advisable to do this with a smaller L, here L=3).
```bash
display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=8,nodeshape=:rect,curves=false))
```
<p align="center"><img src="tutorialFigures/plotTriangularLattice.jpg"></p>

Finally, we perform the embedding to compute the $c_{ii^{\prime}}^{(n)}(i\nu_{m})$ of Eq. (10) in [1], they are provided as vectors containing the prefactors in 
$$c_{ii^{\prime}}^{(n)}(i\nu_{m})=c_{ii^{\prime},0}^{(n)}\delta_{0,m}+(1-\delta_{0,m})\sum_{l=2,4,6,...}c_{ii^{\prime},l}^{(n)}\frac{1}{(2\pi m)^{2l}},$$
c.f. Eq. (17) in [1].
```bash
c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs);
```
Here, the site $i$ is pinned to one of the central basis sites (here the single site 19 in the L=3 example in the figure) and $i^{\prime}$ takes on all other site indices. If this embedding function is used with a SimpleGraph instead where translation symmetry is not assumed, $i^{\prime}$ takes all possible site indices. However, the latter case is less efficient since the embedding with Dyn_HTE_Lattice structures automatically uses lattice symmetries. 

### Equal-time correlators (crosschecks)

As a first crosscheck for Dyn-HTE we reproduce the HTE of the uniform susceptibility $\sum_{i^{\prime}}\left\langle S_{i}^{z}S_{i^{\prime}}^{z}\right\rangle$  found by Elstner et al in [PhysRevLett.71.10 (1993)]. As a first step we analytically sum the $c_{ii^{\prime}}^{(n)}(i\nu_{m})$ over Matsubara frequency to obtain the HTE of the equal-time correlators $\left\langle S_{i}^{z}S_{i^{\prime}}^{z}\right\rangle$, this is done as follows:
```bash
c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat)
```
We now sum over i to obtain the expansion coefficients of the uniform susceptibility in powers of (-x).
```bash
println( [sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max]' )
```
This yields 
```bash
[1/4, 3/8, 3/8, 17/64, 75/512, 441/5120, 8143/122880, 23691/573440, 118351/13762560, -585353/123863040, 46090313/9909043200, 23370989/2076180480, 1154027593/581330534400]
```
which indeed agrees to the result of Elstner et al if their convention for expansion coefficients is taken into account.

We next consider the equal-time correlators in k-space, say at the K-point. We define a vector of inverse temperatures (x_vec) and obtain the Fourier transform to momentum space using the “get_c_k” function. The series expansion in x (instead of -x) is obtained by a simple sign-flip of the even-index entries (note the julia convention that the first element - here $x^{0}$ coefficient - is at index 1). Then the polynomial is obtained as p_x 
```bash
k,k_label = K,"K"
x_vec = collect(0:0.05:5.1)
coeffs_x = flipEvenIndexEntries(get_c_k(k , c_iipEqualTime_mat,hte_lattice))
p_x = Polynomial(coeffs_x)
```
<p align="center"><img src="tutorialFigures/Triangular_EqualTime_GkK.jpg"></p>

The evaluation of the bare series (p_x) is shown in the figure (full green line). It diverges around x=1.5. For a better estimate, we evaluate Padé approximants using, e.g. for [6,6],
```bash
get_pade(p_x,6,6)
```

which provides a rational function that agrees well down to x=5 with the results of the exponential tensor renormalization group (XTRG, geometry YC6x12, D*=1000) by Chen et al in [PhysRevB.99.140404 (2019)] (gray dots). The series in u=\mathrm{tanh}(fx) is obtained as follows from a linear transformation of the vector of expansion coefficients (we pick f=0.2 empirically for good agreement of the u-Padés, blue lines)
```bash
f=0.2 
ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
u_vec = tanh.(f .* x_vec)
p_u = Polynomial(ufromx_mat*coeffs_x)
```
This completes the crosschecking of the frequency-summed Dyn-HTE expansion.

### Static structure factor

We proceed to the study of the static susceptibility $\chi_{\mathbf{k}}\equiv G_{\mathbf{k}}(i\nu_{m}=0)$ at Matsubara index m=0. We obtain its real-space version using the function
```bash
TGiip_Matsubara_xpoly(c_iipDyn_mat,i,1,m)
```

and then compute the spatial Fourier transform by hand (using the cosine due to inversion symmetry and the function for the real-space position of lattice site i
```bash
getSitePosition(hte_lattice.lattice,i)
```
as follows:
```bash
p_x = sum([cos(dot(k,getSitePosition(hte_lattice.lattice,i) 
.- getSitePosition(hte_lattice.lattice,hte_lattice.basis_positions[1]))) 
* get_TGiip_Matsubara_xpoly(c_iipDyn_mat,i,1,m) for i in 1:hte_lattice.lattice.length])
```
For k=K the resulting x-Padés and u-Padés (f=0.25) which are computed from the x-series as above are shown in the top panel of the figure. There we also compare to the bold line diagrammatic Monte Carlo of Kulagin et al [PhysRevB.87.024407(2013)], see dots. Finally we can repeat the above calculation of $\chi_{\mathbf{k}}$ for $\mathbf{k}$ sampled uniformly along a path through the BZ (see bottom panel of figure). This path ($Γ\rightarrow K\rightarrow M\rightarrow Γ$) with Nk+1 $\mathbf{k}$-points and the tick labels at the points defining the polygon is obtained conveniently with:
```bash
path = [Γ,K,M,Γ]
pathticks = ["Γ","K","M","Γ"]
Nk = 200
k_vec,kticks_positioins = create_brillouin_zone_path(path, Nk)
```
<p align="center"><img src="tutorialFigures/Triangular_StaticSF.jpg"></p>





## Short description of scripts

* **Dyn-HTE:** The main code which shows a usage example of Dyn-HTE both for real- and momentum-space Matsubara correlators

* **LatticeGraphs** Generates lattices or lattice balls L by keyword geometry=chain,square,triangle,... (can be expanded by user). The site indices of the center unit-cell of the lattice ball are returned. This code builds on Lattice.jl which is taken from SpinMC.jl

* **Embedding:**  For given L, calculation of embedding factors for graphG and calculation of expansion coefficients of $TG_{ii\prime}(i\nu_m)$ (expansion in powers of -x)

* **ConvenienceFunctions:**  Definition of various functions to help evaluate and plot the results of the Dyn-HSE

* **Dyn-HTE_Tutorial_TriangularLattice.pdf** Tutorial on how to use the Dyn-HTE working through the example of the S=1/2 nearest-neighbor Heisenberg AFM on the triangular lattice.
