push!(LOAD_PATH,"../src/")

using Splittings
using Documenter
using Literate
using Plots # to not capture precompilation output

examples = [
"Vlasov-Ampere"  => "examples/vlasov-ampere",
"Vlasov-Poisson" => "examples/vlasov-poisson",
"Bump On Tail"   => "examples/bump_on_tail",
"Rotation 2D"    => "examples/rotation2d_bsl",
"Vlasov-HMF"     => "examples/vlasov-hmf"
]

# generate examples

for example in examples

    EXAMPLE     = joinpath(@__DIR__, "..",  string(example[2],".jl"))
    DOC_OUTPUT  = joinpath(@__DIR__, "src", "examples")
    NB_OUTPUT   = joinpath(@__DIR__, "src", "notebooks")
   
    Literate.markdown(EXAMPLE, DOC_OUTPUT)
    Literate.notebook(EXAMPLE, NB_OUTPUT, execute=false)

end

makedocs(modules=[Splittings],
         doctest = false,
         format = :html,
         sitename = "Splittings.jl",
         pages = ["Introduction"    => "index.md",
                  "Semi-Lagrangian" => "bsl.md",
		      "Advection functions" => "advections.md",
                  "Examples" => [string(example[2],".md") for example in examples],
                  "User Documentation" => [
                      "How to Contribute" => "contributing.md"],
		          "Contents" => "contents.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pnavaro/Splittings.jl.git",
 )
