module RunExamples

using Test, FaultDetectionTools

@testset "Case Studies" begin
# case studies
include("CS1_1.jl")
include("CS1_2.jl")
include("CS2_1.jl")
include("CS2_2.jl")
end

end
