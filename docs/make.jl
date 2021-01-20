using MeasureTheory
using Documenter

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
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cscherrer/MeasureTheory.jl",
)
