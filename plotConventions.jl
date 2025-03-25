using Plots, LaTeXStrings, ColorSchemes
pyplot()    #using gr plotting backend
Plots.default(titlefont=(7,), legendfontsize=6, background_color_legend=nothing, foreground_color_legend=nothing,
grid=false, guidefont=(7,), tickfont=(7,), framestyle=:box, linewidth=1,markerstrokewidth=0,markersize=3.5,markeralpha=0.9, dpi=200)
aps_width = 325 #in pixels, seems to work with pyplot backend  (243 # width of single-colum plot in points, it assumes internally dpi=100 (so we set this above), but this does not limit the resolution since we later save as pdf
color_vec = [palette(:tab10)[c%10+1] for c in 0:50] #
#color_vec = ["blue","green","red","cyan","brown","purple","magenta","orange","teal","aquamarine2","steelblue1","darkorchid","grey56","olive","blue","green","red","cyan","brown","purple","magenta","orange","teal","aquamarine2","steelblue1","darkorchid","grey56","olive"]
grey_vec = ["black","grey","darkblue","olive","brown","steelblue1"]
thermalCol4_vec = reverse(ColorSchemes.thermal[1:70:end])
marker_vec = [:dot,:cross,:diamond,:dtriangle,:square,:star4,:utriangle]
linestyle_vec = [:solid, :dash, :dot, :dashdot, :solid, :dash, :dot, :dashdot]

plt_empty = plot([0,0],[0,0],label="")

function addABC(plt,label);
    """ adds "label" on the top left of plot axes """
    xl,yl = xlims(plt),ylims(plt)
    annotate!(plt, xl[1]-0.16*(xl[2]-xl[1]), yl[2]-0.02*(yl[2]-yl[1]), (label,8,:left) )
end
