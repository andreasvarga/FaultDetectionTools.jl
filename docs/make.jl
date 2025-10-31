import Pkg; Pkg.add("DescriptorSystems")
using Documenter, DescriptorSystems, FaultDetectionTools, DescriptorSystems
DocMeta.setdocmeta!(FaultDetectionTools, :DocTestSetup, :(using FaultDetectionTools); recursive=true)

makedocs(warnonly = true, 
  modules  = [FaultDetectionTools],
  sitename = "FaultDetectionTools.jl",
  authors  = "Andreas Varga",
  format   = Documenter.HTML(prettyurls = false),
  pages    = [
    "Home"   => "index.md",
    "Tutorials"   => [
      "FDDbasics.md",
      "MDbasics.md",
      "SynthesisParadigms.md"
      ],
    "Fault Detection" => [ 
        "FDIObjects.md",
        "FDIanalysis.md",
        "FDIsynthesis.md",
        "FDIperformance.md"],
    "Model Detection" => [ 
        "MDObjects.md",
        "MDanalysis.md",
        "MDsynthesis.md",
        "MDperformance.md"
        ],
    "Residual Evaluation" => [ 
        "FDDsystems.md"
        ],
     "Utilities" => [
        "FDIutils.md"
     ],
     "Index" => "makeindex.md"
  ]
)

deploydocs(
  repo = "github.com/andreasvarga/FaultDetectionTools.jl.git",
  target = "build",
  devbranch = "master"
)
