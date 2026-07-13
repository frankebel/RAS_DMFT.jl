using Changelog
using Documenter
using Literate
using RAS_DMFT

DocMeta.setdocmeta!(RAS_DMFT, :DocTestSetup, :(using RAS_DMFT); recursive = true)

# generate changelog
Changelog.generate(
    Changelog.Documenter(),
    joinpath(@__DIR__, "../CHANGELOG.md"),
    joinpath(@__DIR__, "src/changelog.md");
    repo = "frankebel/RAS_DMFT.jl",
)

# generate documentation
Literate.markdown("src/tutorial.jl", "src/generated")

makedocs(;
    modules = [RAS_DMFT],
    authors = "Frank Ebel and contributors",
    sitename = "RAS_DMFT.jl",
    format = Documenter.HTML(;
        canonical = "https://frankebel.github.io/RAS_DMFT.jl", edit_link = "main", assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "generated/tutorial.md",
        "API reference" => "api.md",
        "Changelog" => "changelog.md",
    ],
)

# only works for public repos, see <https://github.com/frankebel/RAS_DMFT.jl/settings/pages>
# deploydocs(; repo="github.com/frankebel/RAS_DMFT.jl", devbranch="main")
