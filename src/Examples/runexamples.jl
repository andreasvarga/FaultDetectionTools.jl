module RunExamples

using Test, FaultDetectionTools

@testset "Examples" begin
# solving EFDP
include("Ex5_2.jl")
include("Ex5_3.jl")
include("Ex5_4.jl")
include("Ex5_4c.jl")
# solving AFDP
include("Ex5_5.jl")
include("Ex5_6.jl")
include("Ex5_6c.jl")
include("Ex5_7.jl")
include("Ex5_9.jl")
# solving EFDIP
include("Ex5_10.jl")
include("Ex5_10c.jl")
# solving AFDIP
include("Ex5_11.jl")
include("Ex5_11c.jl")
end

end
