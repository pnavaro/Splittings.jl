import Splittings: UniformMesh, advection!

@testset "BSL advections" begin
  
  p = 5
  nx, nv = 128, 128
  xmin, xmax = -10, 10
  vmin, vmax = -10, 10
  meshx = UniformMesh(xmin, xmax, nx; endpoint=false)
  meshv = UniformMesh(vmin, vmax, nv; endpoint=false)

  f  = zeros(Complex{Float64},(nx,nv))
  f .= exp.(-meshx.x.^2) * transpose(exp.(-meshv.x.^2))
  fᵗ = zeros(Complex{Float64},(nv,nx))

  dt =  0.5

  e = ones(Float64, nx)
  v = ones(Float64, nv)

  advection!(f,  p, meshx, v, nv, dt)
  transpose!(fᵗ, f)
  advection!(fᵗ, p, meshv, e, nx, dt)
  transpose!(f,  fᵗ)
  advection!(f,  p, meshx, -v, nv, dt)
  transpose!(fᵗ, f)
  advection!(fᵗ, p, meshv, -e, nx, dt)
  transpose!(f,  fᵗ)

  f0 =  exp.(-meshx.x.^2) * transpose(exp.(-meshv.x.^2))
  println( maximum( abs.(real(f) -f0)))

  @test real(f) ≈ f0 atol=1e-6

end
