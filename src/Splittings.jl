module Splittings

using Statistics, FFTW, LinearAlgebra

export Mesh, RectMesh2D, UniformMesh
export advection!, advection_v!, advection_x!, advection_y!
export compute_rho, compute_e, interpolatefunction compute_e(mesh

include("meshes.jl")
include("advections.jl")

"""

   Compute charge density
   ρ(x,t) = ∫ f(x,v,t) dv

"""
function compute_rho(meshv, f)
   dv = meshv.dx
   rho = dv * sum(f, dims=2)
   rho .- mean(rho)
end

"""

   Compute Ex using that -ik*Ex = rho

"""
function compute_e(meshx, rho)
   nx = meshx.nx
   k =  2* pi / (meshx.xmax - meshx.xmin)
   modes = zeros(Float64, nx)
   modes .= k * vcat(0:div(nx,2)-1,-div(nx,2):-1)
   modes[1] = 1.0
   rhok = fft(rho)./modes
   real(ifft(-1im*rhok))
end

end # module
