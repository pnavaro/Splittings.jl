# # Vlasov-Poisson
# 
#md # [notebook](@__NBVIEWER_ROOT_URL__notebooks/vlasov-poisson.ipynb),
#
# We consider the dimensionless Vlasov-Poisson equation for one species
# with a neutralizing background.
# 
# ```math
#  \\frac{∂f}{∂t}+ v⋅∇_x f + E(t,x) ⋅ ∇_v f = 0, \\
#  - Δϕ = 1 - ρ, E = - ∇ ϕ \\
#  ρ(t,x)  =  ∫ f(t,x,v)dv.
# ```
# 
#  - [Vlasov Equation - Wikipedia](https://en.wikipedia.org/wiki/Vlasov_equation)
#  - [Landau damping - Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)
# 

using Plots, LinearAlgebra
pyplot()
import Splittings:UniformMesh
import Splittings:@Strang
import Splittings

#-

function push_t!(f, p, meshx, v, nv, dt)
    Splittings.advection!(f, p, meshx, v, nv, dt)
end

#-
#
function push_v!(f, fᵗ, p, meshx, meshv, nrj, dt)
    rho = Splittings.compute_rho(meshv, f)
    e   = Splittings.compute_e(meshx, rho)
    push!(nrj, 0.5*log(sum(e.*e)*meshx.step))
    transpose!(fᵗ, f)
    Splittings.advection!(fᵗ, p, meshv, e, meshx.length, dt)
    transpose!(f, fᵗ)
end

#-

function landau(tf, nt)

  p = 3
  nx, nv = 32, 64
  xmin, xmax = 0.0, 4π
  vmin, vmax = -6., 6.
  meshx = UniformMesh(xmin, xmax, nx; endpoint=false)
  meshv = UniformMesh(vmin, vmax, nv; endpoint=false)
  x = meshx.points
  v = meshv.points
  dx = meshx.step

  ϵ, kx = 0.001, 0.5
  f = zeros(Complex{Float64},(nx,nv))
  f .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) * transpose(exp.(-0.5*v.^2))
  fᵗ = zeros(Complex{Float64},(nv,nx))

  dt = tf / nt

  nrj = Float64[]

  for it in 1:nt
      @Strang( push_t!(f, p, meshx, v, nv, dt),
               push_v!(f, fᵗ, p, meshx, meshv, nrj, dt))
  end

  nrj

end

#-

nt = 600
tf = 60.0
t  = range(0.0, stop=tf, length=nt)
@time nrj = landau(tf, nt)
plot( t, nrj)
plot!(t, -0.1533*t.-5.50)
savefig("landau-plot.png"); nothing # hide

@testset "Vlasov-Poisson" begin  #src
@test length(nrj) > 0            #src
@end                             #src

 
 
#md # ![](landau-plot.png)


