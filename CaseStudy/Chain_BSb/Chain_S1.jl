using JLD2,DelimitedFiles
include("../../plotConventions.jl")
include("../../LatticeGraphs.jl")
include("../../Embedding.jl")
include("../../ConvenienceFunctions.jl") 

### load graph evaluations and prepare lattice  
L = 12
n_max = 1*L
spin_length = 1
hte_graphs = load_dyn_hte_graphs(spin_length,L);

hte_lattice = getLattice(L,"chain");
#display(graphplot(hte_lattice.graph,names=1:nv(hte_lattice.graph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

### compute all correlations in the lattice
c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs)
#c_iiDyn = 1.0 .* c_iipDyn_mat[hte_lattice.basis_positions,1][1]

#########################################################################################
###### Dynamic structure factor (DSF) ###################################################
#########################################################################################

###### dynamical Matsubara correlator (k-space)
k,label=(π,0),"π"
c_kDyn = get_c_k(k,c_iipDyn_mat,hte_lattice)

#label="bulk"
m_vec = get_moments_from_c_kDyn(c_kDyn)
poly_x = Polynomial([0,1],:x)

x_vec_bare = collect(0.0:0.025:1.05)
x_vec = collect(0.0:0.2:6.0)

x0_vec = [1.0,2.0,3.0,4.0,5.0]  # for these x the DSF will be computed

##### plot DSF and related quantities
w_vec = collect(0.0:0.01:3.0)

###### Pade for moments with x-series and u-series
if true
    r_max = 3   # maximal order for moment
    f=0.55 #for k=π
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)
    u_vec = tanh.(f .* x_vec)
    u0_vec = tanh.(f .* x0_vec)
    m0_vec = [Float64[] for _ in x0_vec]

    plt_m = plot([0],[0],xlims=(0,x_vec[end]),label="")
    plot!(plt_m,xlabel=L"x=J/T",ylabel=L"x \cdot m_{\mathbf{k},2r}(x) \, / \, m_{\mathbf{k},2r}(0)",legend=:topleft)
    plot!(plt_m,-x_vec,0*x_vec,color=:grey,label="x bare")
    plot!(plt_m,-x_vec,0*x_vec,color=:grey,linestyle=linestyle_vec[2],label="u Padé [7-r,6-r]")
    plot!(plt_m,-x_vec,0*x_vec,color=:grey,linestyle=linestyle_vec[3],label="u Padé [6-r,5-r]")
    annotate!(plt_m,3,1,text(L"\mathrm{label}="*string(label)*",  f="*string(f),7))
        
    for r in 0:r_max
        println("r=$r")

        xm_norm_r = coeffs(poly_x * (m_vec[1+r]/m_vec[1+r](0)))
        p_x = Polynomial(xm_norm_r)
        p_u = Polynomial(ufromx_mat[1:n_max+2-2*r,1:n_max+2-2*r]*xm_norm_r)
        
        plot!(plt_m,x_vec_bare,Polynomial(xm_norm_r).(x_vec_bare),color=color_vec[r+1],linewidth=0.4,label="r=$r",alpha=0.7) 

        ### x-Padé moments are not well behaved (grow large or negative)
        #plot!(plt_m,x_vec,get_pade(p_x,7-r,6-r).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="",alpha=0.7)
        #plot!(plt_m,x_vec,get_pade(p_x,6-r,5-r).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="",alpha=0.7)

        plot!(plt_m,x_vec,get_pade(p_u,7-r,6-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="",alpha=0.7)
        plot!(plt_m,x_vec,get_pade(p_u,6-r,5-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="",alpha=0.7)

        ### extract moments from u-Padé [7-r,6-r] at x0_vec 
        for x0_pos in eachindex(x0_vec)
            x0 = x0_vec[x0_pos]
            u0 = u0_vec[x0_pos]
            push!(m0_vec[x0_pos],m_vec[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
        end
    end
    display(plt_m)
end


###### δ_r, JS and A for x ∈ x0_vec
ΔHaldane=0.4105
plt_δ=plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_{\mathbf{k},r}",legend=:bottomright)
plt_JS = plot(xlims=(0,w_vec[end]),xlabel=L"\omega/J=w",ylabel=L"J \, S(\mathbf{k}="*label*L",\omega)",legendfontsize=5.0,legend=:topright)
plt_JAo2π = plot(xlabel=L"\omega/J=w",ylabel=L"J \, A(\mathbf{k}="*label*L",\omega)",legend=:topright)
plt_wmax = plot(x0_vec,[ΔHaldane for _ in x0_vec],label="Δ",xlabel="x",ylabel=L"w_{max}",color=:grey)

plot!(plt_JAo2π,[ΔHaldane,ΔHaldane],[0.0,0.3],label="Haldane gap 0.4105",color=:grey)
plot!(plt_JS,[ΔHaldane,ΔHaldane],[0.0,0.3],label="Haldane gap 0.4105",color=:grey)


for x0_pos in eachindex(x0_vec)
    x0 = x0_vec[x0_pos]

    ### plot Dyn-HTE
    δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
    scatter!(plt_δ,r_vec,δ_vec,color=thermalCol13_vec[x0_pos],label="")
    δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,4000,true)
    plot!(plt_δ,r_max+1:6,δ_vec_ext[r_max+2:7],label="",color=thermalCol13_vec[x0_pos])

    JSw_vec = [JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec]

    wmax = w_vec[argmax(JSw_vec)]
    scatter!(plt_wmax,[x0],[wmax],color=thermalCol13_vec[x0_pos],label="")

    plot!(plt_JS,w_vec, JSw_vec,color=thermalCol13_vec[x0_pos],label="")


    plot!(plt_JAo2π,w_vec, JSw_vec .* (1 .- exp.(-x0 .* w_vec)) ,color=thermalCol13_vec[x0_pos],label="")
    #scatter!(plt_JSw0,[1/x0],[JSw_vec[1]],color=thermalCol15_vec[x0_pos],label="")
end

### plot results
xPlots,yPlots=4,1
plt_final = plot(plt_m, plt_δ, plt_JS, plt_wmax, layout=(yPlots,xPlots), size=(aps_width*xPlots,0.6*aps_width*yPlots))
display(plt_final)
savefig(plt_final,"CaseStudy/Chain_BSb/Chain_S1_DSF_k"*label*".png")



##########################################################
###### singlet gap estimation at k=π (from Matsubara correlator, does not work well)
##########################################################
k=(π,0)
k_label="π"
c_kDyn = get_c_k(k,c_iipDyn_mat,hte_lattice)
poly_x =Polynomial([0,1])

function Gkπ_poly(m::Int,c_kDyn)
    if m==0
        poly = Polynomial(flipEvenIndexEntries(c_kDyn[:,1]))
    else
        Δ2l_vec = [1/(2*π*m)^twol for twol in 2:2:18]
        tmp = [sum( c_kDyn[r,2:end] .* Δ2l_vec ) for r in 1:n_max+1]
        #@show tmp
        poly = Polynomial(flipEvenIndexEntries(tmp))
    end
    return poly
end

function gapn(n::Int,x_vec::Vector{Float64},c_kDyn,N,M)
    X=[1/prod([(k+j)*(k-j) for j in 0:n if j != k]) for k in 0:n]

    ### Padé for entire numerator and entire denominator
    num = Polynomial(sum([-k^2*X[1+k]*Gkπ_poly(k,c_kDyn) for k in 0:n]).coeffs[3:end]) #divide out the 1/x^2 from before the sqrt
    den = sum([X[1+k]*Gkπ_poly(k,c_kDyn) for k in 0:n])
    @show num
    @show den 
    return 2*π*1 * real.(Complex.( get_pade(num,N-1,M-1).(x_vec) ./ get_pade(den,N,M).(x_vec) ) .^ 0.5)

    ### Padé for each G
    #num = sum([-k^2*X[1+k]*get_pade(Gkπ_poly(k,c_kDyn),N,M).(x_vec) for k in 0:n])
    #den = sum([X[1+k]*get_pade(Gkπ_poly(k,c_kDyn),N,M).(x_vec) for k in 0:n])
    #return 2*π*1 ./ x_vec .* real.(Complex.( num ./ den ) .^ 0.5)
end

function gapn_test(n::Int,x_vec::Vector{Float64},c_kDyn)
    ### expand the fraction in the gap estimator to order 10, then resum.
    ### does not work well, the argument of the square root becomes negative already for moderate x
    X=[1/prod([(k+j)*(k-j) for j in 0:n if j != k]) for k in 0:n]

    num = -sum([k^2*X[1+k]*Gkπ_poly(k,c_kDyn) for k in 0:n]).coeffs

    den = sum([X[1+k]*Gkπ_poly(k,c_kDyn) for k in 0:n]).coeffs

    @variables x
    tmp = (series(num[3:end],x) // series(den,x))
    @show tmp

    p = [taylor_coeff(tmp,x,m;rationalize=false).val for m in 0:10]
    @show p

    return 2*π*1 * ( get_pade(Polynomial(p),5,5).(x_vec) ) .^ 0.5
end

if true
    x_vec = collect(1.0:0.1:6)
    plt_gap=plot([x_vec[1]],[0],label="",legend=:bottomleft,xlabel=L"J/T=x",ylabel=L"\Delta E_{k=\pi} \; / \; J")

    for n in 1:4     
        plot!(plt_gap,x_vec,gapn(n,x_vec,c_kDyn,6,6),label="n=$n [6,6]")
        #plot!(plt_gap,x_vec,gapn_test(n,x_vec,c_kDyn),label="n=$n")
    end

    plot!(plt_gap,[2.5,x_vec[end]],[0.4105,0.4105],linestyle=:dash,color=:gray,label="DMRG")

    xPlots,yPlots=1,1
    plt_final = plot(plt_gap, layout=(yPlots,xPlots), size=(aps_width*xPlots,0.5*aps_width*yPlots))
    display(plt_final)
    savefig(plt_final,"CaseStudy/Chain_BSb/Chain_S1_HaldaneGap.png")
end