module RunExamples

using Test, FaultDetectionTools

@testset "Examples" begin
# LPV to noise inputs
include("Ex2_4.jl")
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
# solving EMDP
include("Ex6_1.jl")
include("Ex6_1c.jl")
# solving AMDP
include("Ex6_2.jl")
include("Ex6_2c.jl")
# synthesis paradigms: nullspace and least order synthesis
include("Ex7_3.jl")
include("Ex7_3c.jl")
include("Ex7_4.jl")
include("Ex7_4c.jl")
end

end
