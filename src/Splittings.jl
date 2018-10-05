module Splittings

using Statistics, FFTW, LinearAlgebra

export RectMesh1D, RectMesh1D1V
export advection!
export compute_rho, compute_e, interpolate

include("cubic_splines.jl")
include("meshes.jl")
include("fft.jl")
include("bsl.jl")

"""

    compute_rho( mesh, f)

    Compute charge density

    ρ(x,t) = ∫ f(x,v,t) dv

    return ρ - ρ̄ if neutralized=true

"""
function compute_rho(mesh, f, neutralized=true)

   dv  = mesh.dx
   ρ = dv * sum(f, dims=2)
   if (neutralized)
       ρ .- mean(ρ)
   else
       ρ
   end

end

"""

    compute_e( mesh, ρ)

    Compute 1d electric field using that -ik * e = ρ

"""
function compute_e(mesh, ρ)

   nx = mesh.nx
   k =  2π / (mesh.xmax - mesh.xmin)
   modes  = zeros(Float64, nx)
   modes .= k * vcat(0:nx÷2-1,-nx÷2:-1)
   modes[1] = 1.0
   ρ̂ = fft(ρ)./modes
   vec(real(ifft(-1im*ρ̂)))

end

end
