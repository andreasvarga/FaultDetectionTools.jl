using Documenter, FaultDetectionTools, DescriptorSystems
DocMeta.setdocmeta!(FaultDetectionTools, :DocTestSetup, :(using FaultDetectionTools); recursive=true)

makedocs(
  modules  = [FaultDetectionTools],
  sitename = "FaultDetectionTools.jl",
  authors  = "Andreas Varga",
  format   = Documenter.HTML(prettyurls = false),
  pages    = [
     "Home"   => "index.md",
     "Library" => [ 
        "FDIObjects.md",
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
  devbranch = "main"
)
