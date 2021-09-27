module RunExamples

using Test, FaultDetectionTools

@testset "Examples" begin
# solving EFDP
include("Ex5_2.jl")
include("Ex5_3.jl")
include("Ex5_4.jl")
include("Ex5_4c.jl")
# solving EFDIP
include("Ex5_10.jl")
include("Ex5_10c.jl")
end

end
