# Generation of Fig6_2 with CairoMakie for Ex6_1
# The required inputs are: 
#    R    - residual internal form
using Makie, CairoMakie, LaTeXStrings
Makie.inline!(false)

N, M = size(R.sys)

titles = [ "Model $i" for i in 1:M]
ylabelr = [ latexstring("\$r^{($i)}\$") for i in 1:N]

fig2 = Figure(size = (1200, 800))

axs = [Axis(fig2[row, col], yticks = WilkinsonTicks(3)) for row in 1:N, col in 1:M]
y = zeros(101,1,2); tout = zeros(101);
for row in 1:N
    for col in 1:M
        axs[row,col].xgridvisible= false         
        axs[row,col].ygridvisible= false         
        axs[row,col].title = row == 1 ? titles[col] : "" 
        axs[row,col].ylabel = col == 1 ? ylabelr[row] : "" 
        row < N ? hidexdecorations!.(axs[row,col], label = false, grid = false) :
                  axs[row, col].xticks = [0,2.,4]; 
        row == 2 && col == 1 &&  (axs[row, col].yticks = [-1,0])
        y[:,:,:], tout[:], _ = stepresp(R.sys[col,row],4)          
        lines!(axs[row, col], tout, y[:,1,1], color = :blue)
        lines!(axs[row, col], tout, y[:,1,2], color = :red)
        xlims!.(axs[row, col], low = 0, high = 4) 
        row == col && ylims!.(axs[row, col], low = -1, high = 1)
    end
end

Label(fig2[end+1, :], text = "Time (seconds)", font = "TeX Gyre Heros Bold",
                   valign = :top,padding = (0, 0, 5, -10))

fig2                   

# comment out next line to save plot
# save("Fig6_2.pdf", fig2, size = (1200, 800))
