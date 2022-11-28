# Generation of Fig6_1 with CairoMakie for Ex6_1
# The required inputs are: 
#    distinf  - the norms of final residual models indexed by residual numbers and model numbers
using Makie, CairoMakie, LaTeXStrings
Makie.inline!(false)

set_theme!()
n, m = size(distinf)
x = 1:n; y = 1:m; z = distinf
set_theme!(colormap = :Hiroshige)
fig1 = Figure(resolution = (800, 600))
ax3d = Axis3(fig1[1, 1]; aspect = (1, 1, 0.2), perspectiveness = 0.1, elevation = 0.87, azimuth = 3.9, 
             title = "Norms of residual models", 
             xlabel = "Model numbers",
             ylabel = "Residual numbers", zlabel = "",
             )
surface!(ax3d, x, y, z; transparency = true)
wireframe!(ax3d, x, y, z, color=:gray)             
ax3d.xticks = 1:n
ax3d.yticks = 1:m

fig1

# comment out next line to save plot
#save("Fig6_1.pdf", fig1, resolution = (800, 600))



