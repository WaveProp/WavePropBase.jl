using WavePropBase
using Documenter

DocMeta.setdocmeta!(WavePropBase, :DocTestSetup, :(using WavePropBase; using WriteVTK); recursive=true)

makedocs(;
    modules=[WavePropBase],
    repo="https://github.com/WaveProp/WavePropBase.jl/blob/{commit}{path}#{line}",
    sitename="WavePropBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://WaveProp.github.io/WavePropBase.jl",
        assets=String[],
    ),
    pages=[
        "WavePropBase" => "index.md",
        "Utils" => "utils.md",
        "Geometry" => "geometry.md",
        "Interpolation" => "interpolation.md",
        "Integration" => "integration.md",
        "Meshes" => "mesh.md",
        "Trees" => "trees.md",
        "IO" => "io.md",
        "References" => "references.md"
    ],
)

deploydocs(;
    repo="github.com/WaveProp/WavePropBase.jl",
    devbranch="main",
)
