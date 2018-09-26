using FFTW, LinearAlgebra

"""
    advection_v!( f, meshx, meshv, E, dt) 

    Advection in υ
    ∂ f / ∂ t − E(x) ∂ f / ∂ υ  = 0

"""
function advection_v!( fᵀ, meshx::UniformMesh, meshv::UniformMesh, E, dt)

    n = meshv.nx
    L = meshv.xmax - meshv.xmin
    k = 2π/L*[0:n÷2-1;-n÷2:-1]
    ek = exp.(-1im * dt * k * transpose(E))

    fft!(fᵀ, 1)
    fᵀ .= fᵀ .* ek
    ifft!(fᵀ, 1)

end

function advection_x!( f, meshx::UniformMesh, meshv::UniformMesh, dt)

    L = meshx.xmax - meshx.xmin
    m = div(meshx.nx,2)
    k = 2*π/L * [0:1:m-1;-m:1:-1]
    k̃ = 2*π/L * [1;1:1:m-1;-m:1:-1]
    v = meshv.x
    ev = exp.(-1im*dt * k * transpose(v))

    fft!(f,1)
    f .= f .* ev
    Ek  = -1im * meshv.dx * sum(f,dims=2) ./ k̃
    Ek[1] = 0.0
    ifft!(f,1)
    real(ifft(Ek))

end
