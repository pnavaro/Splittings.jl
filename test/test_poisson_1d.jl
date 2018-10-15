using FFTW, LinearAlgebra

import Splittings:UniformMesh, compute_rho, compute_e
import Splittings:Landau, distribution

mesh1 = UniformMesh( 0., 4π,  64, endpoint=false)
mesh2 = UniformMesh(-6., 6., 128, endpoint=false)

ϵ, kx = 0.001, 0.5

f = zeros(Complex{Float64},(mesh1.length,mesh2.length))

f .= distribution( mesh1, mesh2, Landau(ϵ, kx))


@testset "Poisson 1D using FFT" begin

ρ = compute_rho(mesh2, f)
@test real(ρ) ≈ ϵ*cos.(kx*mesh1.points)
e = compute_e(mesh1, ρ)
@test real(e) ≈ ϵ*sin.(kx*mesh1.points)/kx

import Splittings:compute_rho!, compute_e!

ex  = zeros(Complex{Float64},mesh1.length)
rho = similar(ex)

compute_rho!(rho, mesh2, f)
@test real(rho) ≈ ϵ*cos.(kx*mesh1.points)
compute_e!(ex, mesh1, rho)
@test real(ex) ≈ ϵ*sin.(kx*mesh1.points)/kx

end
