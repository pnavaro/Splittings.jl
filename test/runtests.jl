using Test
using IntervalSets
using LinearAlgebra
using Literate

include("geometry.jl")
#include("splinepp.jl")
include("linear_pendulum.jl")
include("test_meshes.jl")
include("test_cubic_splines.jl")
include("test_poisson_1d.jl")
include("test_poisson_2d.jl")
include("test_bsl.jl")
include("test_rotation_2d.jl")

examples = ["vlasov-ampere",
            "vlasov-poisson",
            "bump_on_tail",
            "rotation2d_bsl",
            "vlasov-hmf"]

output = joinpath(@__DIR__, "examples")

for example in examples
    filename = joinpath(@__DIR__, "..", "examples", string(example,".jl"))
    Literate.script(filename, output)
    include(joinpath(@__DIR__, "examples", string(example,".jl")))
end
