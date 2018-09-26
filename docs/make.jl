push!(LOAD_PATH,"../src/")
using Documenter, Splittings
makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/pnavaro/Splittings.jl.git",
    julia = "0.7",
    osname = "osx"
)
