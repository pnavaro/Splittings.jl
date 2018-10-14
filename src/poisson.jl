export compute_rho,  compute_e
export compute_rho!, compute_e!

import Statistics:mean
import FFTW:fft, ifft, fft!, ifft!

"""

    compute_rho( mesh, f)

    Compute charge density

    ρ(x,t) = ∫ f(x,v,t) delta2

    return ρ - ρ̄ if neutralized=true

"""
function compute_rho(mesh2::UniformMesh, f, neutralized=true)

   local delta2  = mesh2.step
   ρ = delta2 * sum(f, dims=2)
   if (neutralized)
       ρ .- mean(ρ)
   else
       ρ
   end

end

"""

    compute_rho!(rho, mesh, f)

    Inplace computation of charge density

"""
function compute_rho!(rho::Vector{Complex{Float64}}, 
		      mesh2::UniformMesh, 
		      f::Array{Complex{Float64},2})

   local delta2 = mesh2.step
   rho .= delta2 * vec(sum(f, dims=2))
   rho .= rho .- mean(rho)

end

"""

    compute_e( mesh, ρ)

    Compute 1d electric field using that -ik * e = ρ

"""
function compute_e(mesh1::UniformMesh, ρ)

   n1 = mesh1.length
   k =  2π / (mesh1.stop - mesh1.start)
   modes  = zeros(Float64, n1)
   modes .= k * vcat(0:n1÷2-1,-n1÷2:-1)
   modes[1] = 1.0
   ρ̂ = fft(ρ)./modes
   vec(real(ifft(-1im*ρ̂)))

end

"""

    compute_e!( e, mesh, ρ)

    Inplace computation of electric field

"""
function compute_e!(e::Vector{Complex{Float64}}, 
		    mesh1::UniformMesh,
		    ρ::Vector{Complex{Float64}})

   n1 = mesh1.length
   k =  2π / (mesh1.stop - mesh1.start)
   modes  = zeros(Float64, n1)
   modes .= k * vcat(0:n1÷2-1,-n1÷2:-1)
   modes[1] = 1.0
   fft!(ρ)
   e .= -1im * ρ ./ modes
   ifft!(e)

end

export RectMesh2D

"""

    mesh = RectMesh2D( mesh1, mesh2)

    Data type that represents a 2D Rectangular Mesh 

"""
struct RectMesh2D

    x1min :: Float64
    x1max :: Float64
    n1   :: Int64
    delta1   :: Float64
    ymin :: Float64
    ymax :: Float64
    n2   :: Int64
    delta2   :: Float64

    function RectMesh2D( mesh1::UniformMesh, mesh2::UniformMesh )

        x1min = mesh1.start    
        x1max = mesh1.stop    
        n1   = mesh1.length
        delta1   = mesh1.step

        ymin = mesh2.start    
        ymax = mesh2.stop    
        n2   = mesh2.length
        delta2   = mesh2.step

        new( x1min, x1max, n1, delta1, ymin, ymax, n2, delta2)

    end

end

export meshgrid
"""

    X, Y = meshgrid( x, y )
    
    Utility function to mimic the numpy function 

"""
function meshgrid( x, y )

   n1, n2 = size(x)[1], size(y)[1]
   repeat(x,1,n2), repeat(y',n1,1)

end


export poisson!
"""

   poisson!(ρ, mesh, ex, ey)

   Solve the equation Δ Φ = - ρ

   ex = ∂ Φ / ∂ x
   ey = ∂ Φ / ∂ y

   Be careful the ρ array is destroyed

"""
function poisson!( ρ::Array{Complex{Float64},2}, 
		   mesh::RectMesh2D, 
		   ex::Array{Complex{Float64},2}, 
		   ey::Array{Complex{Float64},2} )

    kx0 =  2π / (mesh.x1max - mesh.x1min)
    ky0 =  2π / (mesh.ymax - mesh.ymin)

    fft!(ρ,[1,2])
    
    kx = kx0 * vcat(0:mesh.n1÷2-1,-mesh.n1÷2:-1)
    kx[1] = 1.0
    ky = ky0 * vcat(0:mesh.n2÷2-1,-mesh.n2÷2:-1)
    kx[1] = 1.0

    for i = 1:mesh.n1
       kx2 = kx[i]*kx[i]
       for j =  1:mesh.n2÷2+1
          k2 = kx2 +ky[j]*ky[j]
          ex[i,j] = -1im * kx[i]/k2 * ρ[i,j]
	  ey[i,j] = -1im * ky[j]/k2 * ρ[i,j]
       end 
       for j = mesh.n2÷2+2:mesh.n2            
          k2 = kx2 +ky[j]*ky[j]
          ex[i,j] = -1im * kx[i]/k2 * ρ[i,j]
          ey[i,j] = -1im * ky[j]/k2 * ρ[i,j]
       end 
    end

    ifft!(ex,[1,2])
    ifft!(ey,[1,2])

end
