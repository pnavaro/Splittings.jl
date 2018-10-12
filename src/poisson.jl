export compute_rho,  compute_e
export compute_rho!, compute_e!

import Statistics:mean
import FFTW:fft, ifft, fft!, ifft!

"""

    compute_rho( mesh, f)

    Compute charge density

    ρ(x,t) = ∫ f(x,v,t) dv

    return ρ - ρ̄ if neutralized=true

"""
function compute_rho(mesh, f, neutralized=true)

   local dv  = mesh.dx
   ρ = dv * sum(f, dims=2)
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
		      mesh, f::Array{Complex{Float64},2})

   local dv = mesh.dx
   rho .= dv * vec(sum(f, dims=2))
   rho .= rho .- mean(rho)

end

"""

    compute_e( mesh, ρ)

    Compute 1d electric field using that -ik * e = ρ

"""
function compute_e(mesh, ρ)

   nx = mesh.nx
   k =  2π / (mesh.xmax - mesh.xmin)
   modes  = zeros(Float64, nx)
   modes .= k * vcat(0:nx÷2-1,-nx÷2:-1)
   modes[1] = 1.0
   ρ̂ = fft(ρ)./modes
   vec(real(ifft(-1im*ρ̂)))

end

"""

    compute_e!( e, mesh, ρ)

    Inplace computation of electric field

"""
function compute_e!(e::Vector{Complex{Float64}}, mesh, ρ::Vector{Complex{Float64}})

   nx = mesh.nx
   k =  2π / (mesh.xmax - mesh.xmin)
   modes  = zeros(Float64, nx)
   modes .= k * vcat(0:nx÷2-1,-nx÷2:-1)
   modes[1] = 1.0
   fft!(ρ)
   e .= -1im * ρ ./ modes
   ifft!(e)

end

export RectMesh2D

"""

    mesh = RectMesh2D( meshx, meshy)

    Data type that represents a 2D Rectangular Mesh 

"""
struct RectMesh2D

    xmin :: Float64
    xmax :: Float64
    nx   :: Int64
    dx   :: Float64
    ymin :: Float64
    ymax :: Float64
    ny   :: Int64
    dy   :: Float64

    function RectMesh2D( meshx::UniformMesh, meshy::UniformMesh )

        xmin = meshx.xmin    
        xmax = meshx.xmax    
        nx   = meshx.nx
        dx   = meshx.dx

        ymin = meshy.xmin    
        ymax = meshy.xmax    
        ny   = meshy.nx
        dy   = meshy.dx

        new( xmin, xmax, nx, dx, ymin, ymax, ny, dy)

    end

end

export meshgrid
"""

    X, Y = meshgrid( x, y )
    
    Utility function to mimic the numpy function 

"""
function meshgrid( x, y )

   nx, ny = size(x)[1], size(y)[1]
   repeat(x,1,ny), repeat(y',nx,1)

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

    kx0 =  2π / (mesh.xmax - mesh.xmin)
    ky0 =  2π / (mesh.ymax - mesh.ymin)

    fft!(ρ,[1,2])
    
    kx = kx0 * vcat(0:mesh.nx÷2-1,-mesh.nx÷2:-1)
    kx[1] = 1.0
    ky = ky0 * vcat(0:mesh.ny÷2-1,-mesh.ny÷2:-1)
    kx[1] = 1.0

    for i = 1:mesh.nx
       kx2 = kx[i]*kx[i]
       for j =  1:mesh.ny÷2+1
          k2 = kx2 +ky[j]*ky[j]
          ex[i,j] = -1im * kx[i]/k2 * ρ[i,j]
	  ey[i,j] = -1im * ky[j]/k2 * ρ[i,j]
       end 
       for j = mesh.ny÷2+2:mesh.ny            
          k2 = kx2 +ky[j]*ky[j]
          ex[i,j] = -1im * kx[i]/k2 * ρ[i,j]
          ey[i,j] = -1im * ky[j]/k2 * ρ[i,j]
       end 
    end

    ifft!(ex,[1,2])
    ifft!(ey,[1,2])

end
