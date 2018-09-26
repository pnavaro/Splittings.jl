push!(LOAD_PATH,"../src/")
using Documenter, Splittings

makedocs(modules=[Splittings],
         doctest = false,
         format = :html,
         sitename = "Splittings",
         pages = ["Introduction" => "index.md",
                  "User Documentation" => [
                    "Semi-Lagrangian" => "bsl.md",
                    "How to Contribute" => "contributing.md"]])


deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/pnavaro/Splittings.jl.git",
    julia = "0.7",
    osname = "osx"
)
