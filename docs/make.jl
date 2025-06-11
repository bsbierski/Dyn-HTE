using DynHTE
using Documenter

DocMeta.setdocmeta!(DynHTE, :DocTestSetup, :(using DynHTE); recursive=true)

makedocs(;
    modules=[DynHTE],
    authors="BenediktSchneiderLMU <Schneider.benedikt@physik.uni-muenchen.de> and contributors",
    sitename="DynHTE.jl",
    format=Documenter.HTML(;
        canonical="https://bsbierski.github.io/Dyn-HTE.jl",
        edit_link="package",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Adding New Lattices" => "Lattices.md",
        "Convenience Functions" => "convenience_functions.md",
         "Embedding of Graphs" => "Embedding.md"
    ],
    checkdocs = :exports,
)

deploydocs(;
    repo="github.com/bsbierski/DynHTE.jl",
    devbranch="package",
)

