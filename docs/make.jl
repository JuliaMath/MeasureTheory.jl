using Documenter
using MeasureBase
using MeasureTheory

DocMeta.setdocmeta!(MeasureBase, :DocTestSetup, :(using MeasureBase); recursive=true)
DocMeta.setdocmeta!(MeasureTheory, :DocTestSetup, :(using MeasureTheory); recursive=true)

pages = [
    "Home" => "index.md",
    "Tutorials" => [
        "Adding a New Measure" => "adding.md",
        "Affine transformations" => "affine.md",
    ],
    "MeasureBase API" => "api_measurebase.md",
    "MeasureTheory API" => "api_measuretheory.md",
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
