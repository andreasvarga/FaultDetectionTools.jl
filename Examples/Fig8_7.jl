# Generation of Fig8_7 with CairoMakie for CS2_1.jl
# The required inputs are: 
#    y    - the collection of step responses
#    tout - the time samples  
using Makie, CairoMakie, LaTeXStrings
Makie.inline!(false)

N = length(y)
title_u = [ L"From: $f_1$"  L"From: $f_2$"  L"From: $u_1$"  L"From: $u_2$" L"From: $u_3$" L"From: $d_1$"  L"From: $d_2$"]
ylabel_r = [ latexstring("To: \$r_1\$")  latexstring("To: \$r_2\$") ]
ns, pp, mm = size(y[1])
fig3 = Figure(resolution = (800, 500))

axs = [Axis(fig3[row, col]) for row in 1:pp, col in 1:mm]
inputs = [mu+md .+ (1:mf); 1:mu+md]
for row in 1:pp
    for col in 1:mm
        #axs[row,col].xgridvisible= false         
        #axs[row,col].ygridvisible= false  
        axs[row,col].yminorgridvisible = false       
        axs[row,col].xminorgridvisible = false       
        axs[row,col].title = row == 1 ? title_u[col] : "" 
        axs[row,col].ylabel = col == 1 ? ylabel_r[row] : "" 
        col > 1 ? hideydecorations!(axs[row,col], label = false, grid = false) :
                  axs[row, col].yticks = [-1,0,1,2] 
        row < 2 ? hidexdecorations!.(axs[row,col], label = false, grid = false) :
                  axs[row, col].xticks = [0,5,10]; 
        #band!(axs[row, col], tout, yl[:,row,col],  yu[:,row,col], color = :skyblue)
        for k = 1:N
            lines!(axs[row, col], tout, y[k][:,row,inputs[col]], color = k == nom ? (:black, 3) : (:blue, 3))
        end
        xlims!.(axs[row, col], low = 0, high = 10) 
        ylims!.(axs[row, col], low = -1, high = 2)
    end
end

Label(fig3[end+1, :], text = "Time (seconds)", font = "TeX Gyre Heros Bold",
                   valign = :top,padding = (0, 0, 5, -10))

fig3

# comment out next line to save plot
# save("Fig8_7.pdf", fig3, resolution = (800, 500))
