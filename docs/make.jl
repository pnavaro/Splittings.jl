push!(LOAD_PATH,"../src/")

using Splittings
using Documenter
using Literate
using Plots # to not capture precompilation output

examples = [
"Vlasov-Ampere"  => "examples/vlasov-ampere.md",
"Vlasov-Poisson" => "examples/vlasov-poisson.md",
"Bump On Tail"   => "examples/bump_on_tail.md",
"Rotation 2D"    => "examples/rotation2d_bsl.md",
"Vlasov-HMF"     => "examples/vlasov-hmf.md"
]

# generate examples

for example in examples

    EXAMPLE     = joinpath(@__DIR__, "..",  string(example[2][1:end-3],".jl"))
    DOC_OUTPUT  = joinpath(@__DIR__, "src", "examples")
    NB_OUTPUT   = joinpath(@__DIR__, "..",  "notebooks")
   
    Literate.markdown(EXAMPLE, DOC_OUTPUT)
    Literate.notebook(EXAMPLE, NB_OUTPUT)

end

makedocs(modules=[Splittings],
         doctest = true,
         format = :html,
         sitename = "Splittings.jl",
         pages = ["Introduction"    => "index.md",
                  "Semi-Lagrangian" => "bsl.md",
		  "Advection functions" => "advections.md"
                  "Examples" => examples,
                  "User Documentation" => [
                    "How to Contribute" => "contributing.md"],
		  "Contents" => "contents.md"
		  ])

deploydocs(
    repo   = "github.com/pnavaro/Splittings.jl.git",
    julia  = "1.0",
    osname = "osx"
 )
