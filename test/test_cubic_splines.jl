import Splittings:compute_interpolants, interpolate

@testset "CubicSpline interpolation" begin

    function interpolation_test(nx = 128)
    
        xmin, xmax = 0.0, 1.0
        x = collect(range(xmin,stop=xmax,length=nx))
        y = sin.(2π*x)
        x_new = zeros(Float64,nx)
        y_new = zeros(Float64,nx)
        coeffs = compute_interpolants(nx, y) 
        for (i, xi) in enumerate(x)
            x_new[i] = xi - 0.1
            x_new[i] = xmin + mod(x_new[i]-xmin,xmax-xmin) 
            y_new[i] = interpolate(coeffs, nx, xmin, xmax, x_new[i])
        end
        maximum(abs.(sin.(2π*(x.-0.1)) - y_new))
    
    end

    @test ≈(interpolation_test(), 0.0, atol=1e-7)

end


import LinearAlgebra: transpose
import Splittings: UniformMesh, advection!

@testset " BSL advections with cubic splines " begin
  
  nx, nv = 128, 128
  xmin, xmax = -5, 10
  vmin, vmax = -5, 10
  meshx = UniformMesh(xmin, xmax, nx)
  meshv = UniformMesh(vmin, vmax, nv)

  f  = zeros(Float64,(nx,nv))
  f .= exp.(-meshx.points.^2) * transpose(exp.(-meshv.points.^2))
  fᵗ = zeros(Float64,(nv,nx))

  dt =  0.5

  e = ones(Float64, nx)
  v = ones(Float64, nv)

  advection!(f, meshx,  v, dt)
  advection!(f, meshv,  e, dt)
  advection!(f, meshx, -v, dt)
  advection!(f, meshv, -e, dt)

  f0 =  exp.(-meshx.points.^2) * transpose(exp.(-meshv.points.^2))
  println( maximum( abs.(f .- f0)))

  @test f ≈ f0 atol=1e-3

end
