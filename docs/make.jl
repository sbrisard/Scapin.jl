using Dates
using Documenter
using Scapin

DocMeta.setdocmeta!(Scapin, :DocTestSetup, :(using Scapin); recursive = true)

use_KaTeX = false

# Use MathJax syntax to define macros:
# - macro names should not be preceded by "\\"
macros = Dict(
    "cellindices" => "\\mathcal{P}",
    "conj" => "\\operatorname{conj}",
    "D" => "\\mathrm d",
    "dbldot" => "\\mathbin{\\mathord{:}}",
    "dft" => "\\operatorname{DFT}",
    "element" => "\\mathrm{e}",
    "fftfreq" => "\\operatorname{Z}",
    "I" => "\\mathrm i",
    "integers" => "\\mathbb Z",
    "naturals" => "\\mathbb N",
    "PI" => "\\mathrm{\\pi}",
    "reals" => "\\mathbb R",
    "sinc" => "\\operatorname{sinc}",
    "symotimes" => "\\stackrel{\\mathrm{s}}{\\otimes}",
    "strain" => "\\boldsymbol{\\varepsilon}",
    "strains" => "\\mathcal E",
    "stress" => "\\boldsymbol{\\sigma}",
    "stresses" => "\\mathcal S",
    #"tens" => "\\@ifstar{\\boldsymbol}{\\mathbf}",
    "tens" => "\\boldsymbol",
    "tensors" => "\\mathcal T",
    "tr" => "\\operatorname{tr}",
    "tuple" => "\\mathsf",
    "vec" => "\\boldsymbol",
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
                #"packages" => ["base", "ams", "autoload", "bm"],
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

open(joinpath(@__DIR__, "src", "assets", "custom.sty"), "w") do io
    write(io, "% WARNING: Automatically generated -- Do NOT edit by hand!\n")
    write(io, "% Script: $(@__FILE__)\n")
    write(io, "% Date: $(now())\n\n")
    for (name, definition) in
        macros write(io, "\\providecommand{\\$name}{$definition}\n")
    end
end

makedocs(;
    modules = [Scapin],
    authors = "Sébastien Brisard <sbrisard@users.noreply.github.com> and contributors",
    repo = "https://github.com/sbrisard/Scapin.jl/blob/{commit}{path}#{line}",
    sitename = "Scapin.jl",
    format = format,
    pages = [
        "Home" => "index.md",
        "Convolution operators" => "convolution_operators.md",
        "Tutorials" => [],
        "Library" => "api.md",
        "Design" => "design.md",
        "Appendix" => "appendix.md",
    ],
)

deploydocs(; repo = "github.com/sbrisard/Scapin.jl", versions = nothing)
