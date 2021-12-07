module RunExamples

using Test, FaultDetectionTools

@testset "Examples and Case Studies" begin
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
# solving EMMP
include("Ex5_12.jl")
include("Ex5_12c.jl")
include("Ex5_13.jl")
include("Ex5_13a.jl")
include("Ex5_13c.jl")
# solving AMMP
include("Ex5_14.jl")
include("Ex5_15.jl")
include("Ex5_16.jl")
include("Ex5_16c.jl")
include("Ex5_17.jl")
end

end
