push!(LOAD_PATH,"../src/")
using Documenter, Splittings

makedocs(modules=[Splittings],
         doctest = true,
         format = :html,
         sitename = "Splittings",
         pages = ["Introduction"    => "index.md",
                  "Semi-Lagrangian" => "bsl.md",
                  "Examples" => [
                      "Vlasov-Ampere"  => "examples/vlasov-ampere.md",
                      "Vlasov-Poisson" => "examples/vlasov-poisson.md",
                      "Bump On Tail"   => "examples/bump_on_tail.md",
                      "Rotation 2D"    => "examples/rotation2d_bsl.md"],
                  "User Documentation" => [
                    "How to Contribute" => "contributing.md"]])


deploydocs(
    repo   = "github.com/pnavaro/Splittings.jl.git",
    julia  = "0.7",
    osname = "osx"
 )
