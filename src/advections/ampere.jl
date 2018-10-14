import FFTW: fft!, ifft!
import LinearAlgebra: transpose

export Ampere

struct Ampere end

"""
    advection!( fᵀ, mesh1, mesh2, E, dt, type, axis ) 

    if axis == 1 Advection in x and compute electric field

    ∂ f / ∂ t − υ ∂ f / ∂ x  = 0

    ∂E / ∂t = −J = ∫ fυ dυ

    if axis == 2 Advection in υ

    ∂ f / ∂ t − E(x) ∂ f / ∂ υ  = 0

"""
function advection!( f     :: Array{Complex{Float64},2}, 
                     fᵀ    :: Array{Complex{Float64},2},
		     mesh1 :: UniformMesh, 
		     mesh2 :: UniformMesh, 
		     e     :: Vector{Complex{Float64}}, 
		     dt    :: Float64,
	             type  :: Ampere, 
		     axis  :: Int64 )

    @assert ( axis == 1 || axis == 2 )

    if ( axis == 2 ) 

        n = mesh2.length ÷ 2
        L = mesh2.stop - mesh2.start
        k = 2π/L * [0:1:n-1;-n:1:-1]
        ek = exp.(-1im * dt * k * transpose(e))

        transpose!(fᵀ, f)
        fft!(fᵀ, 1)
        fᵀ .= fᵀ .* ek
        ifft!(fᵀ, 1)
        transpose!(f, fᵀ)

    else 

        L = mesh1.stop - mesh1.start
        m = mesh1.length ÷ 2
        k = 2π/L * [0:1:m-1;-m:1:-1]
        v = mesh2.points
        ev = exp.(-1im * dt * k * transpose(v))

        fft!(f,1)
        f   .= f .* ev
        k[1] = 1.0
        e   .= -1im * mesh2.step * vec(sum(f,dims=2)) ./ k
        e[1] = 0im
        ifft!(f,1)
        ifft!(e)

     end 

end
