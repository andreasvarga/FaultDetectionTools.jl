# Generation of Fig8_4 with CairoMakie for CS1_2.jl
# The required inputs are: 
#    y    - the collection of step responses
#    tout - the time samples  
using Makie, CairoMakie, LaTeXStrings
Makie.inline!(false)

N = length(y)
title_u = reshape([ latexstring("From: \$f_$i\$")  for i in 1:mf],1,mf)
ylabel_r = reshape([ latexstring("To: \$r_$i\$")  for i in 1:mf],1,mf)
ns, pp, mm = size(y[1])
fig = Figure(resolution = (800, 600))

axs = [Axis(fig[row, col]) for row in 1:pp, col in 1:mm]

for row in 1:pp
    for col in 1:mm
        axs[row,col].xgridvisible= false         
        axs[row,col].ygridvisible= false  
        axs[row,col].yminorgridvisible = false       
        axs[row,col].xminorgridvisible = false       
        axs[row,col].title = row == 1 ? title_u[col] : "" 
        axs[row,col].ylabel = col == 1 ? ylabel_r[row] : "" 
        col > 1 ? hideydecorations!(axs[row,col], label = false, grid = false) :
                  axs[row, col].yticks = [0, 1]
        row < 8 ? hidexdecorations!.(axs[row,col], label = false, grid = false) :
                  axs[row, col].xticks = [0,1,2]; 
        #band!(axs[row, col], tout, yl[:,row,col],  yu[:,row,col], color = :skyblue)
        for k = 1:N
            lines!(axs[row, col], tout, y[k][:,row,col], color = (:blue, 3))
        end
        xlims!.(axs[row, col], low = 0, high = 2) 
        ylims!.(axs[row, col], low = -0.1, high = 1.1)
    end
end

Label(fig[end+1, :], text = "Time (seconds)", font = "TeX Gyre Heros Bold",
                   valign = :top,padding = (0, 0, 5, -10))

fig

# comment out next line to save plot
# save("Fig8_4.pdf", fig, resolution = (800, 600))
