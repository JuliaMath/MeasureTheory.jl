using MeasureTheory
using Documenter


pages = [
    "Introduction" => "intro.md"
    "Home" => "index.md"
    "Adding a New Measure" => "adding.md"
]

makedocs(;
    modules=[MeasureTheory],
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
