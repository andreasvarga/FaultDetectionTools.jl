module Runtests

using Test, FaultDetectionTools

@testset "Test FaultDetectionTools" begin
# test FDI tools
include("test_fdiutils.jl")
# test FDI synthesis functions
include("test_efdsyn.jl")
end

end
