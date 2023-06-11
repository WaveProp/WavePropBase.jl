using WavePropBase
using Documenter

# go over all .md files in the src folder and replace ```julia by ```@example.
# This is a hack so that e.g. vscode can properly highlight the code in the
# markdown files, and documenter can still run the examples.
dir  = joinpath(WavePropBase.PROJECT_ROOT,"docs","src")
for file in readdir(dir)
    fname,fext = splitext(file)
    if fext == ".md"
        old_content = read(joinpath(dir,file), String)
        new_content = replace(old_content, r"```julia" => "```@example")
        write(joinpath(dir,fname*".md"), new_content)
    end
end

DocMeta.setdocmeta!(WavePropBase,
                    :DocTestSetup,
                    :(using WavePropBase);
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
                "Examples" => ["helmholtz_scattering_2d.md"],
                # "Utils" => "utils.md",
                # "Geometry" => "geometry.md",
                # "Interpolation" => "interpolation.md",
                # "Integration" => "integration.md",
                # "Meshes" => "mesh.md",
                # "Trees" => "trees.md",
                # "IO" => "io.md",
                "References" => "references.md"
                ])

for file in readdir(dir)
    fname,fext = splitext(file)
    if fext == ".jmd"
        old_content = read(joinpath(dir,file), String)
        new_content = replace(old_content, "```@example" => "```julia")
        write(joinpath(dir,fname*".md"), new_content)
    end
end

deploydocs(; repo="github.com/WaveProp/WavePropBase.jl", devbranch="main")
