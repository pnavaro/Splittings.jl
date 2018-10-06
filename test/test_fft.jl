using FFTW, LinearAlgebra

import Splittings:RectMesh1D, compute_rho, compute_e

meshx = RectMesh1D( 0., 2π,  64, endpoint=false)
meshv = RectMesh1D(-6., 6., 128, endpoint=false)

ϵ, kx = 0.001, 0.5

f = zeros(Float64,(meshx.nx,meshv.nx))

f .= (1.0.+ϵ*cos.(kx*meshx.x))/sqrt(2π) .* transpose(exp.(-0.5*meshv.x.^2))

ρ = compute_rho(meshv, f)
@test ρ ≈ ϵ*cos.(kx*mesh.x)
e = compute_e(meshx, ρ)
@test e ≈ ϵ*sin.(kx*mesh.x)/kx
