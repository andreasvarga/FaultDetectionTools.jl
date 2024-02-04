# Generation of Fig8_3 with CairoMakie for CS1_1.jl
# The required inputs are: 
#    y    - the collection of step responses
#    tout - the time samples  
using Makie, CairoMakie, LaTeXStrings
Makie.inline!(false)

N = length(y)
title_u = reshape([ latexstring("From: \$f_$i\$")  for i in 1:mf],1,mf)
ylabel_r = reshape([ latexstring("To: \$r_$i\$")  for i in 1:6],1,6)
#yticks = [[-50,0,50],[-50,0,50],[-50,0,50],[-50,0,50],[-10,0,10],[0,1,2]]
yticks = [[-50,0,50],[-50,0,50],[-50,0,50],[-50,0,50],[-10,0,10],["   0","   1","   2"]]
yhighs = [50,50,50,51,10,2]
ylows = [-50,-51,-51,-50,-10,-0.2]
ns, p, m = size(y[1])
fig2 = Figure(size = (800, 600))

axs = [Axis(fig2[row, col]) for row in 1:p, col in 1:m]

for row in 1:p
    for col in 1:m
        axs[row,col].xgridvisible= false         
        #axs[row,col].ygridvisible= false  
        axs[row,col].yminorgridvisible = false       
        axs[row,col].title = row == 1 ? title_u[col] : "" 
        axs[row,col].ylabel = col == 1 ? ylabel_r[row] : "" 
        col > 1 ? hideydecorations!(axs[row,col], label = false, grid = false) :
                  (row < 6 ? axs[row, col].yticks = yticks[row] : axs[row, col].yticks = ([0,1,2],yticks[row])) 
        row < 6 ? hidexdecorations!.(axs[row,col], label = false, grid = false) :
                  axs[row, col].xticks = [0,5.,10]; 
        #band!(axs[row, col], tout, yl[:,row,col],  yu[:,row,col], color = :skyblue)
        for k = 1:N
            lines!(axs[row, col], tout, y[k][:,row,col], color = (:blue, 3))
        end
        xlims!.(axs[row, col], low = 0, high = 10) 
        ylims!.(axs[row, col], low = ylows[row], high = yhighs[row])
    end
end

Label(fig2[end+1, :], text = "Time (seconds)", font = "TeX Gyre Heros Bold",
                   valign = :top,padding = (0, 0, 5, -10))

fig2

# comment out next line to save plot
#save("Fig8_3.pdf", fig2, size = (800, 600))
