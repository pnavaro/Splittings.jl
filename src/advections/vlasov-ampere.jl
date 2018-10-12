import FFTW: fft!, ifft!
import LinearAlgebra: transpose

export advection_v!
"""
    advection_v!( fᵀ, meshx, meshv, E, dt) 

    Advection in υ

    ∂ f / ∂ t − E(x) ∂ f / ∂ υ  = 0

    Apply this advection on transposed `f`

"""
function advection_v!( fᵀ    :: Array{Complex{Float64},2}, 
		       meshx :: UniformMesh, 
		       meshv :: UniformMesh, 
		       e     :: Vector{Complex{Float64}}, 
		       dt    :: Float64)

    n = meshv.nx
    L = meshv.xmax - meshv.xmin
    k = 2π/L*[0:n÷2-1;-n÷2:-1]
    ek = exp.(-1im * dt * k * transpose(e))

    fft!(fᵀ, 1)
    fᵀ .= fᵀ .* ek
    ifft!(fᵀ, 1)

end

export advection_x!
"""
    advection_x!( f, meshx, meshv, dt) 

    Advection in x and compute electric field

    ∂ f / ∂ t − υ ∂ f / ∂ x  = 0

    ∂E / ∂t = −J = ∫ fυ dυ

"""
function advection_x!( f     :: Array{Complex{Float64},2}, 
		       meshx :: UniformMesh, 
		       meshv :: UniformMesh, 
		       e     :: Vector{Complex{Float64}}, 
		       dt    :: Float64)

    L = meshx.xmax - meshx.xmin
    m = meshx.nx ÷ 2
    k = 2π/L * [0:1:m-1;-m:1:-1]
    v = meshv.x
    ev = exp.(-1im*dt * k * transpose(v))

    fft!(f,1)
    f   .= f .* ev
    k[1] = 1.0
    e   .= -1im * meshv.dx * vec(sum(f,dims=2)) ./ k
    e[1] = 0.0im
    ifft!(f,1)
    ifft!(e)

end
