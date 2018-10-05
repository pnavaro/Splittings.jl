using FFTW, LinearAlgebra


meshx = Splittings.RectMesh1D(xmin, xmax, nx, endpoint=false)
meshv = Splittings.RectMesh1D(vmin, vmax, nv, endpoint=false)

x = meshx.x
v = meshv.x
ϵ, kx = 0.001, 0.5

f = zeros(Complex{Float64},(nx,nv))
fᵀ= zeros(Complex{Float64},(nv,nx))

f .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))
transpose!(fᵀ,f)

ρ = Splittings.compute_rho(meshv, f)
e = Splittings.compute_e(meshx, ρ)

dt = 0.1

advection!(fᵀ, meshx, meshv, e,  0.5dt)
transpose!(f,fᵀ)
e = Splittings.advection!( f, meshx, meshv, dt)
transpose!(fᵀ,f)
advection!(fᵀ, meshx, meshv, e,  0.5dt)
