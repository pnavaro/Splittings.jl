module Splittings

using Statistics, FFTW, LinearAlgebra

export advection!
export interpolate

include("meshes.jl")
include("domains.jl")
include("fft.jl")
include("bsl.jl")
include("cubic_splines.jl")
include("operator_splitting.jl")
include("poisson.jl")

end
