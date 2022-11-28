# Generation of Fig6_3 with CairoMakie for Ex6_2
# The required inputs are: 
#    R    - residual internal form
using CairoMakie, LaTeXStrings
using LinearAlgebra

N, M = size(R.sys)

# generate input signals for Ex. 6.2
d = diagm([ 1, 1, 0.01, 0.01, 0.01, 0.01, 0.03, 0.03])
t = Vector(0:0.01:10);  ns = length(t);
usin = sin.(2*t).+1.5
tau = 2*pi; usquare = +(rem.(t,tau) .>= tau/2) .+ 0.3
u = [ usquare usin (rand(ns,mw).-0.5)]*d;  

titles = [ "Model $i" for i in 1:M]
ylabelr = [ latexstring("\$\\theta_$i\$") for i in 1:N]
#ylabelr = [ L"\theta_$i" for i in 1:N]
# ind = [2, 4, 5, 9]
# ylabelr[ind] = [ L"\theta_$i\newline" for i in ind]

fig = Figure(;font = "CMU Serif", fontsize=14, resolution = (1200, 800))

axs = [Axis(fig[row, col], yticks = WilkinsonTicks(3)) for row in 1:N, col in 1:M]
nr = size(R.sys[1,1],1)
r = zeros(ns,nr); tout = zeros(ns); 
α = 0.9; β = 0.1; γ = 10;
s = rtf('s'); filt = dss(rtf(1/(s+ γ)))
for row in 1:N
    for col in 1:M
        axs[row,col].xgridvisible= false         
        axs[row,col].ygridvisible= false         
        axs[row,col].title = row == 1 ? titles[col] : "" 
        axs[row,col].ylabel = col == 1 ? ylabelr[row] : "" 
        row < N ? hidexdecorations!.(axs[row,col], label = false, grid = false) :
                  axs[row, col].xticks = [0,5.,10]; 
        r[:,:], tout[:], _ = timeresp(R.sys[col,row],u,t)          
        θ = α*sqrt.(r[:,1].^2 .+ r[:,2].^2) +
            β*sqrt.(timeresp(filt,r[:,1].^2 .+ r[:,2].^2,t)[1])
        lines!(axs[row, col], tout, θ, color = :blue)
        xlims!.(axs[row, col], low = 0, high = 10) 
        ylims!.(axs[row, col], low = 0, high = 8)
        col == 1 ? axs[row, col].yticks = [0,4,8] : hideydecorations!.(axs[row,col], label = false, grid = false) 
    end
end

Label(fig[end+1, :], text = "Time (seconds)", font = "TeX Gyre Heros Bold",
                   valign = :top,padding = (0, 0, 5, -10))

colgap!(fig.layout, 20)
rowgap!(fig.layout, 15)


fig                   

# comment out next line to save plot
#save("Fig6_3.pdf", fig, resolution = (1200, 800))
