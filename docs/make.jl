using WavePropBase
using Documenter
using Weave

# write .md files from .jmd files
weave_files = ["index.jmd","utils.jmd"]
for fname in weave_files
    dir  = joinpath(WavePropBase.PROJECT_ROOT,"docs","src")
    weave(joinpath(dir, fname), doctype="github")
end

DocMeta.setdocmeta!(WavePropBase,
                    :DocTestSetup,
                    :(using WavePropBase; using WriteVTK);
                    recursive=true)

makedocs(;
         modules=[WavePropBase],
         repo="https://github.com/WaveProp/WavePropBase.jl/blob/{commit}{path}#{line}",
         sitename="WavePropBase.jl",
         format=Documenter.HTML(;
                                prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://WaveProp.github.io/WavePropBase.jl",
                                assets=String[]),
         pages=["WavePropBase" => "index.md",
                # "Utils" => "utils.md",
                # "Geometry" => "geometry.md",
                # "Interpolation" => "interpolation.md",
                # "Integration" => "integration.md",
                # "Meshes" => "mesh.md",
                # "Trees" => "trees.md",
                # "IO" => "io.md",
                # "References" => "references.md"
                ])

deploydocs(; repo="github.com/WaveProp/WavePropBase.jl", devbranch="main")
