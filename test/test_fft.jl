using FFTW, LinearAlgebra

import Splittings:RectMesh1D, compute_rho, compute_e, advection!

meshx = RectMesh1D( 0., 2π,  64, endpoint=false)
meshv = RectMesh1D(-6., 6., 128, endpoint=false)

ϵ, kx = 0.001, 0.5

f = zeros(Complex{Float64},(meshx.nx,meshv.nx))
fᵀ= zeros(Complex{Float64},(meshv.nx,meshx.nx))

f .= (1.0.+ϵ*cos.(kx*meshx.x))/sqrt(2π) .* transpose(exp.(-0.5*meshv.x.^2))
transpose!(fᵀ,f)

ρ = compute_rho(meshv, f)
e = compute_e(meshx, ρ)

dt = 0.1

advection!(fᵀ, meshx, meshv, e,  0.5dt)
transpose!(f,fᵀ)
e = advection!( f, meshx, meshv, dt)
transpose!(fᵀ,f)
advection!(fᵀ, meshx, meshv, e,  0.5dt)
