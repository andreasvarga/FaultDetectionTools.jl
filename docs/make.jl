import Pkg; Pkg.add("DescriptorSystems")
using Documenter, DescriptorSystems, FaultDetectionTools
DocMeta.setdocmeta!(FaultDetectionTools, :DocTestSetup, :(using FaultDetectionTools); recursive=true)

makedocs(
  modules  = [FaultDetectionTools],
  sitename = "FaultDetectionTools.jl",
  authors  = "Andreas Varga",
  format   = Documenter.HTML(prettyurls = false),
  pages    = [
    "Home"   => "index.md",
    "Tutorials"   => [
      "FDDbasics.md",
      "SynthesisParadigms.md"
      ],
    "Library" => [ 
      "FDIObjects.md",
        "FDIanalysis.md",
        "FDIperformance.md",
        "FDIsynthesis.md"
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
