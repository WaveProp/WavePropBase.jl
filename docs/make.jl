using WavePropBase
using Documenter

DocMeta.setdocmeta!(WavePropBase, :DocTestSetup, :(using WavePropBase); recursive=true)

makedocs(;
    modules=[WavePropBase],
    authors="Luiz M. Faria <maltezfaria@gmail.com> and contributors",
    repo="https://github.com/WaveProp/WavePropBase.jl/blob/{commit}{path}#{line}",
    sitename="WavePropBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://WaveProp.github.io/WavePropBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/WaveProp/WavePropBase.jl",
    devbranch="main",
)