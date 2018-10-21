module Splittings

export UniformMesh
export advection!

include("meshes.jl")
include("domains.jl")
include("geometry.jl")
include("operator_splitting.jl")
include("poisson.jl")
include("advections/ampere.jl")
include("advections/bspline.jl")
include("advections/cubic_splines.jl")
include("advections/splinepp.jl")
include("advections/splinenn.jl")
include("initializers/landau.jl")
include("initializers/hmf.jl")
include("initializers/bump_on_tail.jl")

end
