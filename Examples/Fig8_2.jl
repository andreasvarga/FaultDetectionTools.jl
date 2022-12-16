# Generation of Fig8_2 with CairoMakie for CS1_1.jl
# The required inputs are: 
#    y    - the nominal step responses
#    tout - the time samples  
using Makie, CairoMakie, LaTeXStrings
Makie.inline!(false)

title_u = reshape([ latexstring("From: \$f_$i\$")  for i in 1:mf],1,mf)
ylabel_r = reshape([ latexstring("To: \$r_$i\$")  for i in 1:6],1,6)
yticks = [[-20,0,20],[-20,0,20],[-20,0,20],[0,10,20],[0,1,2],[0,1,2]]
yhighs = [20,20,20,20,2,2]
ylows = [-20,-20,-20,-1,-0.1,-0.1]
ns, p, m = size(y)
f = Figure(resolution = (800, 600))

axs = [Axis(f[row, col]) for row in 1:p, col in 1:m]

for row in 1:p
    for col in 1:m
        axs[row,col].xgridvisible= false         
        #axs[row,col].ygridvisible= false  
        axs[row,col].yminorgridvisible = false       
        axs[row,col].title = row == 1 ? title_u[col] : "" 
        axs[row,col].ylabel = col == 1 ? ylabel_r[row] : "" 
        col > 1 ? hideydecorations!(axs[row,col], label = false, grid = false) :
                  axs[row, col].yticks = yticks[row]
        row < 6 ? hidexdecorations!.(axs[row,col], label = false, grid = false) :
                  axs[row, col].xticks = [0,5.,10]; 
        #band!(axs[row, col], tout, yl[:,row,col],  yu[:,row,col], color = :skyblue)
        lines!(axs[row, col], tout, y[:,row,col], color = (:blue, 3))
        xlims!.(axs[row, col], low = 0, high = 10) 
        ylims!.(axs[row, col], low = ylows[row], high = yhighs[row])
    end
end

Label(f[end+1, :], text = "Time (seconds)", font = "TeX Gyre Heros Bold",
                   valign = :top,padding = (0, 0, 5, -10))

f

# comment out next line to save plot
#save("Fig8_2.pdf", f, resolution = (800, 600))