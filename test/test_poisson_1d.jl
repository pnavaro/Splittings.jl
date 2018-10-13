using FFTW, LinearAlgebra

import Splittings:UniformMesh, compute_rho, compute_e

meshx = UniformMesh( 0., 4π,  64, endpoint=false)
meshv = UniformMesh(-6., 6., 128, endpoint=false)

ϵ, kx = 0.001, 0.5

f = zeros(Complex{Float64},(meshx.length,meshv.length))

f .= (1.0.+ϵ*cos.(kx*meshx.points))/sqrt(2π) .* transpose(exp.(-0.5*meshv.points.^2))

@testset "Poisson 1D using FFT" begin

ρ = compute_rho(meshv, f)
@test real(ρ) ≈ ϵ*cos.(kx*meshx.points)
e = compute_e(meshx, ρ)
@test real(e) ≈ ϵ*sin.(kx*meshx.points)/kx

import Splittings:compute_rho!, compute_e!

ex  = zeros(Complex{Float64},meshx.length)
rho = similar(ex)


compute_rho!(rho, meshv, f)
@test real(rho) ≈ ϵ*cos.(kx*meshx.points)
compute_e!(ex, meshx, rho)
@test real(ex) ≈ ϵ*sin.(kx*meshx.points)/kx

end
