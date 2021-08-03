using Scapin
using Documenter

DocMeta.setdocmeta!(Scapin, :DocTestSetup, :(using Scapin); recursive = true)

use_KaTeX = true

# Use MathJax syntax to define macros:
# - macro names should not be preceded by "\\"
macros = Dict(
    "cellindices" => "\\mathcal P",
    "conj" => "\\operatorname{conj}",
    "D" => "\\mathrm d",
    "dbldot" => "\\mathbin{\\mathord{:}}",
    "dft" => "\\operatorname{DFT}",
    "fftfreq" => "\\operatorname{Z}",
    "I" => "\\mathrm i",
    "integers" => "\\mathbb Z",
    "naturals" => "\\mathbb N",
    "PI" => "\\mathrm{\\pi}",
    "reals" => "\\mathbb R",
    "sinc" => "\\operatorname{sinc}",
    "strains" => "\\mathcal E",
    "stresses" => "\\mathcal S",
    "sym" => "\\operatorname{\\textbf{\\textsf{sym}}}",
    "tens" => "\\bm",
    "tensors" => "\\mathcal T",
    "tuple" => "\\mathsf",
    "vec" => "\\bm",
)

if use_KaTeX
    macros2 = Dict()
    for (key, val) ∈ pairs(macros)
        macros2["\\"*key] = val
    end
    macros = macros2
end

mathengine =
    use_KaTeX ?
    KaTeX(
        Dict{Symbol,Any}(
            :delimiters => Dict{Any,Any}[
                Dict(:left => "\$", :right => "\$", display => false),
                Dict(:left => "\$\$", :right => "\$\$", display => true),
                Dict(:left => "\\[", :right => "\\]", display => true),
            ],
            :macros => macros,
        ),
    ) :
    MathJax3(
        Dict{Symbol,Any}(
            :tex => Dict{String,Any}(
                "macros" => macros,
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
            "theory/appendix.md",
            "theory/bibliography.md"
        ],
        "Library" => "api.md",
    ],
)

deploydocs(; repo = "github.com/sbrisard/Scapin.jl", versions = nothing)
