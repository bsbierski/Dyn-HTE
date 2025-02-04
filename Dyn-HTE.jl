using JLD2
#using RobustPade


import Pkg 
Pkg.activate(@__DIR__) #activates the environment in the folder of the current file

include("Embedding.jl")
include("GraphGeneration.jl")
include("LatticeGraphs.jl")
include("ConvenienceFunctions.jl") 
#specify max order
max_order = 10

#LOAD FILES 
#-------------------------------------------------------------------------------------

#generate list of graphs
graphs_vec = [load_object("GraphFiles/graphs_"*string(nn)*".jld2") for nn in 0:max_order];
gG_vec = getGraphsG(graphs_vec);
## identify which gG have the same underlying simple-graph structure. Precalculate Symmetry factors. 
gG_vec_unique = give_unique_gG_vec(gG_vec);

#create vector of all lower order dictionaries
C_Dict_vec = Vector{Vector{Vector{Rational{Int64}}}}(undef,max_order+1) 
   
#load dictionaries of all lower orders C_Dict_vec 
for ord = 0:max_order
    C_Dict_vec[ord+1]  = load_object("GraphFiles/GraphG_Lists/C_"*string(ord)*".jld2")
end 
#-----------------------------------

#1. Define lattice ball for embedding (it is enough for embedding of max_order graphs to have ball radius L=max_order)
L = 4
lattice,LatGraph,center_sites = getLattice_Ball(L,"honeycomb");
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.1,fontsize=7,nodeshape=:rect,curves=false))

#2.Compute all correlations in the lattice
Correlators = compute_lattice_correlations(LatGraph,lattice,center_sites,max_order,gG_vec_unique,C_Dict_vec);

fourier_transform((1.0*pi,1pi),Correlators,lattice,center_sites)[2,1]
#3. Fourier Transform



### Compute A 2D Brillouin zone cut: 
N = 100
kx = range(-8pi,8pi,length=N)
ky = range(-8pi,8pi,length=N) #[0.] #for chains
kmat = [(y,x,x) for x in kx, y in ky ]
structurefactor =  brillouin_zone_cut(kmat,Correlators,lattice,center_sites);

### Evaluate the correlators at a frequency and plot the 2D Brillouin zone cut
iomega = 0
JoverT = -10.5
padetype = [4,4]
evaluate(x) = eval_correlator_LR_continuous_pad_Mats(x,iomega, JoverT, padetype); #define evaluation function
struc = real.(evaluate.(structurefactor));
p = Plots.heatmap(kx,ky,struc,clims=(0,1.))

############# Brillouin zone path
#1. Define a high symmetry path through the brillouin zone

#---chain
path = [(0,0),(2pi,0)]
pathticks = ["Γ","Γ"]

#---triangular
path = [(0,0),(0,2pi/(sqrt(3))),(2pi/(3),2pi/(sqrt(3))),(0.,0.)]
pathticks = ["Γ","M","K","Γ"]

#---triangular Sherman 1
path = [(0,2pi/(sqrt(3))),(0.001,0.001),(4pi/(3),0),(1.0*pi,pi/sqrt(3)),(2pi/(3),0),(pi/2,pi/(2*sqrt(3)))]
pathticks = ["M","Γ","K","M","Y1","Y"]

#---triangular Sherman 2
path = [(0.1,0),(4pi/(3),0),(1.0*pi,pi/sqrt(3)),(4pi/(3),0),(0.1,0)]
pathticks = ["Γ","K","M","K","Γ"]

#---square
path = [(0,0),(pi,0),(pi,pi),(0,0)]
pathticks = ["Γ","X","M","Γ"]

#---square path Sherman
path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
pathticks = ["K","X","M","K","Γ","X"]


#Generate the path 
Nk = 100
kvec,kticks_positioins = create_brillouin_zone_path(path, Nk)
BrillPath = brillouin_zone_path(kvec,Correlators,lattice,center_sites);

kticks_positioins
kvec[63]
BrillPath

####Plot the susceptibility along the high symmetry path for various temperatures 
#= Tvec = [2.0,1.0,0.5,0.375]
padetype = [4,4]

p= plot();
for T in Tvec
evaluate(x) = eval_correlator_LR_continuous_pad_Mats(x,0,-1/T,padetype);
cut = 1/T* real.(evaluate.(BrillPath));
plot!(p,1:length(kvec),cut, label = "T = $T", xticks=(kticks_positioins,pathticks)#= ,yrange = [0,2],xrange = [900,2800] =#)
end
p

BrillPath[1]

get_intDiffApprox(p::Polynomial,Tvec,4,2,2)


####Plot the Spectrum along the high symmetry path for a set temperature  
ωvec = -6:0.1:6
JoverT = 0.5
padetype = [4,4]
spectrum = transpose([imag(eval_correlator_LR_continuous_pad(x,ω,JoverT,padetype)) for x in BrillPath, ω in ωvec]);
heatmap(1:length(kvec),ωvec,spectrum, xticks=(kticks_positioins,pathticks),colormap = :RdBu )

 =#

###### S(k,w) heatmap}
using CairoMakie

x = 2.0
#k_vec = [(k,0.0) for k in 0.01:(2*π-0.02)/14:(2*π-0.01)]
w_vec = collect(0.01:0.0314/2:5.0)
JSkw_mat = get_JSkw_mat_finitex("pade",x,kvec,w_vec,0.02,0,3,200,false,Correlators,lattice,center_sites)





fig = Figure(fontsize=25,resolution = (900,500));
ax=Axis(fig[1,1],limits=(0,Nk+1,0,5),ylabel=L"\omega/J=w",title="Square Lattice: x="*string(x),titlesize=25,xlabelsize=25,ylabelsize=25);
hm=CairoMakie.heatmap!(ax,[k for k in 1:Nk+1],w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,1.0),highclip=:white);
ax.xticks = (kticks_positioins,pathticks)
CairoMakie.Colorbar(fig[:, end+1], hm,size=40, label = L"J S(k,w)")
resize_to_layout!(fig);
CairoMakie.plot!(ax,QMC_axis*(Nk+0.5),1*QMCdata, color = :pink, alpha = 0.45,label = L"T=0 \text{QMC-Gap}")
CairoMakie.plot!(ax,[k for k in 1:Nk+1],2.4*disp, color = :orange, alpha = 0.45,label = L"T=0 \text{LSWT-Dispersion}")
axislegend(ax)
display(fig)

save("SquareX"*string(x)*".pdf",fig)


#plot(JSkw_mat[40,:])



display(fig)
kvec
disp = [sqrt(1-1/4*(cos(k[1])+cos(k[2]))^2)  for k in kvec]


QMCdata =[2.404673091700819,
2.4075874976713854,
2.3834208504098355,
2.3640636136849604,
2.346919825819672,
2.3117723081222055,
2.280859957153502,
2.257086263505961,
2.2144507963859907,
2.1779716025521605,
2.1242679419709383,
2.124071465163934,
2.1307643525251185,
2.117994921596578,
2.131095511014344,
2.144305728795438,
2.1603083376024586,
2.1662671984777515,
2.1848995456165357,
2.19652462675644,
2.1954238684960794,
2.166061365632318,
2.1479106510804766,
2.082221905271982,
2.0149220496460503,
1.930248518799002,
1.8284507528510332,
1.7014071409479683,
1.5619598627940126,
1.3810148931498825,
1.2245351504284647,
1.0375547224292099,
0.8438285907228016,
0.6404096032041728,
0.42766576405152223,
0.216478145957526,
0.30102866523844973,
0.6033846404619969,
0.8992974238875875,
1.182358129191505,
1.4354857488822648,
1.6774796828427716,
1.9004153665238446,
2.07495795839273,
2.2014833271236953,
2.3238228856184797,
2.3836087847469702,
2.407652989940387,
2.377330069392697,
2.3041752049180326,
2.2018070417713917,
2.051779390146115,
1.88175006985842,
1.6772832060357672,
1.4359722628805616,
1.1762673481743657,
0.8989075889530551,
0.6030790098733232,
0.30718493852458995,
0.23173784463487346,
0.44989571135831374,
0.6542939642324885,
0.8607910883941874,
1.0537968051415794,
1.2416941249068552,
1.410336717585693,
1.5678484720242336,
1.718509165627227,
1.8389836182004466,
1.9537403109408407,
2.0442254738729506,
2.103167996194379,
2.164426709784836,
2.204209052985948,
2.222939841920374,
2.2219878650102456,
2.2155641649590163,
2.210593301741803,
2.180401053864168,
2.185150124732715,
2.1570733313817327,
2.1572105532786883,
2.144082991803278]

QMC_axis=[0.0011862269433817046,
0.014651505760147045,
0.029399192083270985,
0.04430717934468976,
0.05761215722316027,
0.07075683416333597,
0.08390151110351166,
0.09864919742663561,
0.1137174856263492,
0.12846517194947316,
0.1400068395067006,
0.15491482676811935,
0.17126552247419155,
0.18216598627823968,
0.19341376878192662,
0.2033396493451823,
0.212639194648086,
0.2225618227285357,
0.23325319835219926,
0.2427597409536837,
0.2534511165773473,
0.26295765917883174,
0.27337722016799537,
0.28395708209545384,
0.29453694402291236,
0.304956505012076,
0.31425395943317586,
0.32387201573086544,
0.3344518776583239,
0.34503173958578237,
0.35448949494517706,
0.3650693568726356,
0.37308440378737684,
0.3849466732211939,
0.3944044285805886,
0.4048239895697523,
0.4291897321905658,
0.4439374185136897,
0.4572423963921602,
0.47215038365357903,
0.48529506059375466,
0.5003633487934682,
0.5138286276102335,
0.5303396242546006,
0.5436446021330711,
0.5567892790732468,
0.5716972663346656,
0.5866052535960843,
0.6016735417957979,
0.6148182187359736,
0.6278025947378545,
0.6428708829375681,
0.6577788701989868,
0.6712441490157521,
0.686152136277171,
0.7010601235385897,
0.7142048004787654,
0.7276700792955307,
0.7425780665569495,
0.7667835082394682,
0.7775236711052215,
0.7869814264646162,
0.7975612883920746,
0.8071793446897642,
0.8177592066172227,
0.8283390685446812,
0.8377968239040758,
0.8483766858315344,
0.857834441190929,
0.8682540021800927,
0.8795017846837797,
0.8889251898421111,
0.899379101032338,
0.9086498386303887,
0.9189090986812576,
0.9297427704276828,
0.939868446363306,
0.9488853741423899,
0.9589843332549639,
0.9696757088786274,
0.9791822514801118,
0.9896018124692756,
0.9980977621989013]

QMC_axis[1]