using Test
using IntervalSets
using LinearAlgebra
using Literate

include("linear_pendulum.jl")
include("test_meshes.jl")
include("test_domains.jl")
include("test_cubic_splines.jl")
include("test_poisson_1d.jl")
include("test_poisson_2d.jl")
include("test_bsl.jl")

examples = ["vlasov-ampere",
            "vlasov-poisson",
            "bump_on_tail",
            "rotation2d_bsl",
            "vlasov-hmf"]

output = joinpath(@__DIR__, "examples")
for example in examples
    Literate.script(joinpath(@__DIR__, '..', 'examples', string(example,".jl"), output))
end

for example in examples
    include(joinpath(@__DIR__, "examples", string(example,".jl")))
end
