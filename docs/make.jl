using Documenter
using MeasureTheory

DocMeta.setdocmeta!(MeasureBase, :DocTestSetup, :(using MeasureBase); recursive=true)
DocMeta.setdocmeta!(MeasureTheory, :DocTestSetup, :(using MeasureTheory); recursive=true)

pages = [
    "Home" => "index.md",
    "Tutorials" => [
        "Adding a new measure" => "adding.md",
        "Affine transformations" => "affine.md",
    ],
    "API Reference" => [
        "MeasureBase" => "api_measurebase.md",
        "MeasureTheory" => "api_measuretheory.md",
        "Index" => "api_index.md",
    ],
]

makedocs(;
    modules=[MeasureBase, MeasureTheory],
    authors="Chad Scherrer <chad.scherrer@gmail.com> and contributors",
    repo="https://github.com/cscherrer/MeasureTheory.jl/blob/{commit}{path}#L{line}",
    sitename="MeasureTheory.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cscherrer.github.io/MeasureTheory.jl",
        assets=String[],
    ),
    pages=pages
)

deploydocs(;
    repo="github.com/cscherrer/MeasureTheory.jl",
)
