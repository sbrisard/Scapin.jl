using Scapin
using Documenter

DocMeta.setdocmeta!(Scapin, :DocTestSetup, :(using Scapin); recursive = true)

mathengine = MathJax3(
    Dict{Symbol,Any}(
        :tex => Dict{String,Any}(
            "macros" => Dict(
                "cellindices" => "\\mathcal P",
                "conj" => "\\operatorname{conj}",
                "D" => "\\mathrm d",
                "dbldot" => "\\mathbin{\\mathord{:}}",
                "dft" => "\\operatorname{DFT}",
                "fftfreq" => "\\operatorname{Z}",
                "I" => "\\mathrm i",
                "integers" => "\\mathbb Z",
                "naturals" => "\\mathbb N",
                "PI" => "\\symup{\\pi}",
                "reals" => "\\mathbb R",
                "strains" => "\\mathcal E",
                "stresses" => "\\mathcal S",
                "sym" => "\\operatorname{\\symbfsf{sym}}",
                "tens" => ["\\symbfsf{#1}", 1],
                "tensors" => "\\mathcal T",
                "tuple" => ["\\mathsf{#1}", 1],
                "vec" => ["\\symbf{#1}", 1],
            ),
            #"packages" => ["base", "ams", "autoload"],
            "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
            "tags" => "ams",
        ),
        :options => Dict(
            "ignoreHtmlClass" => "tex2jax_ignore",
            "processHtmlClass" => "tex2jax_process",
        ),
    ),
    true,
)

format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://sbrisard.github.io/Scapin.jl",
    assets = String[],
    mathengine = mathengine,
)

# format = Documenter.LaTeX(platform = "none")

makedocs(;
    modules = [Scapin],
    authors = "Sébastien Brisard <sbrisard@users.noreply.github.com> and contributors",
    repo = "https://github.com/sbrisard/Scapin.jl/blob/{commit}{path}#{line}",
    sitename = "Scapin.jl",
    format = format,
    pages = [
        "Home" => "index.md",
        "Theory" => [
            "theory/nomenclature.md",
            "theory/continuous_green_operators.md",
            "theory/discrete_green_operators.md",
            "theory/appendix.md"
        ],
        "Library" => "api.md",
    ],
)

deploydocs(; repo = "github.com/sbrisard/Scapin.jl", versions = nothing)
