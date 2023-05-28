module Runtests

using Test, FaultDetectionTools

@testset "Test FaultDetectionTools" begin
# test FDI tools
# include("test_fdiutils.jl")
# # test FDI analysis tools
# include("test_fdigenspec.jl")
# # test FDI synthesis functions
# include("test_efdsyn.jl")
# include("test_efdisyn.jl")
include("test_afdsyn.jl")
# include("test_afdisyn.jl")
# include("test_emmsyn.jl")
# include("test_ammsyn.jl")
# test model detection functions
include("test_emdsyn.jl")
# include("test_amdsyn.jl")
end

end
