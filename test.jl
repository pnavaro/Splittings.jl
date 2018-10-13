using Literate

examples = ["vlasov-ampere",
            "vlasov-poisson",
            "bump_on_tail",
            "rotation2d_bsl",
            "vlasov-hmf"]

output = joinpath(@__DIR__, "examples")

for example in examples
    filename = joinpath(@__DIR__, "..", "examples", string(example,".jl"))
    println(output)
    println(filename)
end

