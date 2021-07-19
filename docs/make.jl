using Scapin
using Documenter

DocMeta.setdocmeta!(Scapin, :DocTestSetup, :(using Scapin); recursive=true)

makedocs(;
    modules=[Scapin],
    authors="Sébastien Brisard <sbrisard@users.noreply.github.com> and contributors",
    repo="https://github.com/sbrisard/Scapin.jl/blob/{commit}{path}#{line}",
    sitename="Scapin.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
