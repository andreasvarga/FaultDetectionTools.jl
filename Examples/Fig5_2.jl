# Generation of Fig5_2 with CairoMakie for Ex5_16
# The required inputs are: 
#    y    - the nominal step responses
#    yl   - the lower bounds of step responses
#    yu   - the upper bounds of step responses 
#    tout - the time samples  
using Makie, CairoMakie, LaTeXStrings
Makie.inline!(false)

title_u = [ latexstring("From: \$f_1\$")  latexstring("From: \$f_2\$")  latexstring("From: \$u_1\$")  latexstring("From: \$u_2\$")]
ylabel_r = [ latexstring("To: \$r_1\$")  latexstring("To: \$r_2\$") ]

ns, p, m = size(y)
fig = Figure(resolution = (800, 600))

axs = [Axis(fig[row, col]) for row in 1:p, col in 1:m]

for row in 1:p
    for col in 1:m
        axs[row,col].xgridvisible= false         
        axs[row,col].ygridvisible= false         
        axs[row,col].title = row == 1 ? title_u[col] : "" 
        axs[row,col].ylabel = col == 1 ? ylabel_r[row] : "" 
        col > 1 ? hideydecorations!(axs[row,col], label = false, grid = false) :
                  axs[row, col].yticks = [0,0.5,1,1.5]
        row < 2 ? hidexdecorations!.(axs[row,col], label = false, grid = false) :
                  axs[row, col].xticks = [0,5.,10]; 
        band!(axs[row, col], tout, yl[:,row,col],  yu[:,row,col], color = :skyblue)
        lines!(axs[row, col], tout, y[:,row,col], color = (:black, 2))
        xlims!.(axs[row, col], low = 0, high = 10) 
        ylims!.(axs[row, col], low = -0.28, high = 1.5)
    end
end

Label(fig[0, :], text = "Step Responses",  fontsize = 20, 
               font = "TeX Gyre Heros Bold", valign = :bottom, 
               padding = (0, 0, -10, 0))
Label(fig[end+1, :], text = "Time (seconds)", font = "TeX Gyre Heros Bold",
                   valign = :top,padding = (0, 0, 5, -10))
Label(fig[0:end, 0], text = "Residuals", font = "TeX Gyre Heros Bold", rotation = pi/2, 
                   valign = :center, padding = (0, -20, 0, 0))
axs[1,1].yticks = [0,0.5,1,1.5]
axs[2,1].yticks = [0,0.5,1,1.5]

fig

# comment out next line to save plot
# save("Fig5_2.pdf", fig, resolution = (800, 600))
