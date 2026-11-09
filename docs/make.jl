using CLASSmodel
using Documenter

DocMeta.setdocmeta!(CLASSmodel, :DocTestSetup, :(using CLASSmodel); recursive=true)

makedocs(;
    modules=[CLASSmodel],
    checkdocs=:exports,
    authors="Thomas Struys",
    sitename="CLASSmodel.jl",
    format=Documenter.HTML(;
        canonical="https://EnvDroneSense.github.io/CLASSmodel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EnvDroneSense/CLASSmodel.jl",
    devbranch="main",
)
