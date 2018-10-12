module Splittings

include("meshes.jl")
include("domains.jl")
include("operator_splitting.jl")
include("poisson.jl")
include("advections/vlasov-ampere.jl")
include("advections/bspline.jl")
include("advections/cubic_splines.jl")

end
