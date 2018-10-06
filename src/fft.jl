using FFTW, LinearAlgebra

"""
    advection!( f, meshx, meshv, E, dt) 

    Advection in υ
    ∂ f / ∂ t − E(x) ∂ f / ∂ υ  = 0

"""
function advection!( fᵀ, meshx::UniformMesh, meshv::UniformMesh, e, dt::Float64)

    n = meshv.nx
    L = meshv.xmax - meshv.xmin
    k = 2π/L*[0:n÷2-1;-n÷2:-1]
    ek = exp.(-1im * dt * k * transpose(e))

    fft!(fᵀ, 1)
    fᵀ .= fᵀ .* ek
    ifft!(fᵀ, 1)

end

"""
    advection!( f, meshx, meshv, dt) 

    Advection in x
    ∂ f / ∂ t − υ ∂ f / ∂ x  = 0

"""
function advection!( f, meshx::UniformMesh, meshv::UniformMesh, dt::Float64)

    L = meshx.xmax - meshx.xmin
    m = div(meshx.nx,2)
    k = 2π/L * [0:1:m-1;-m:1:-1]
    k̃ = 2π/L * [1;1:1:m-1;-m:1:-1]
    v = meshv.x
    ev = exp.(-1im*dt * k * transpose(v))

    fft!(f,1)
    f .= f .* ev
    Ek  = -1im * meshv.dx * sum(f,dims=2) ./ k̃
    Ek[1] = 0.0
    ifft!(f,1)
    real(ifft(Ek))

end
