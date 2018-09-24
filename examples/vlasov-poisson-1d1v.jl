
using Plots, LinearAlgebra, Splittings


"""
    landau(tf, nt)

 Compute Landau damping by solving Vlasov-Poisson system in 1D1V space.
 
 ## Semi-Lagrangian method

 Let us consider an abstract scalar advection equation of the form
 
 ``\\frac{∂f}{∂t}+ a(x, t) ⋅ ∇f = 0.``
 
 The characteristic curves associated to this equation are the solutions of 
 the ordinary differential equations
 
 ``\\frac{dX}{dt} = a(X(t), t)``

 We shall denote by ``X(t, x, s)`` the unique solution of this equation 
 associated to the initial condition ``X(s) = x``.

 The classical semi-Lagrangian method is based on a backtracking of 
 characteristics. Two steps are needed to update the distribution function 
 ``f^{n+1}`` at ``t^{n+1}`` from its value ``f^n`` at time ``t^n`` :

 1. For each grid point ``x_i`` compute ``X(t^n; x_i, t^{n+1})`` the value 
    of the characteristic at ``t^n`` which takes the value ``x_i`` at 
    ``t^{n+1}``.
 2. As the distribution solution of first equation verifies
    ``f^{n+1}(x_i) = f^n(X(t^n; x_i, t^{n+1})),``
    we obtain the desired value of ``f^{n+1}(x_i)`` by computing 
    ``f^n(X(t^n;x_i,t^{n+1})`` by interpolation as ``X(t^n; x_i, t^{n+1})`` 
    is in general not a grid point.

 *[Eric Sonnendrücker - Numerical methods for the Vlasov equations](http://www-m16.ma.tum.de/foswiki/pub/M16/Allgemeines/NumMethVlasov/Num-Meth-Vlasov-Notes.pdf)*


 Vlasov-Poisson equation
 -----------------------

 We consider the dimensionless Vlasov-Poisson equation for one species
 with a neutralizing background.

 ``
 \\frac{∂f}{∂t}+ v⋅∇_x f + E(t,x) ⋅ ∇_v f = 0, \\
 - Δϕ = 1 - ρ, E = - ∇ ϕ \\
 ρ(t,x)  =  ∫ f(t,x,v)dv.
 ``

 - [Vlasov Equation - Wikipedia](https://en.wikipedia.org/wiki/Vlasov_equation)
 - [Landau damping - Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)

# Examples
```julia
nt = 600
tf = 60.0
t  = range(0.0, stop=tf, length=nt)
@time nrj = landau(tf, nt)
plot( t, nrj)
plot!(t, -0.1533*t.-5.50)
```

"""
function landau(tf, nt)

  # Set grid
  p = 3
  nx, nv = 64, 128
  xmin, xmax = 0.0, 4π
  vmin, vmax = -6., 6.
  meshx = UniformMesh(xmin, xmax, nx)
  meshv = UniformMesh(vmin, vmax, nv)
  x = meshx.x
  v = meshv.x
  dx = meshx.dx

  # Initialize distribution function for Landau damping.

  eps, kx = 0.001, 0.5
  f = zeros(Complex{Float64},(nx,nv))
  f .= (1.0.+eps*cos.(kx*x))/sqrt(2.0*pi) * transpose(exp.(-0.5*v.^2))
  fᵗ = zeros(Complex{Float64},(nv,nx))

  # Set time domain
  dt = tf / nt

  # Run simulation
  nrj = Float64[]

  for it in 1:nt
     advection!(f, p, meshx, v, nv, 0.5*dt)
     rho = compute_rho(meshv, f)
     e   = compute_e(meshx, rho)
     transpose!(fᵗ, f)
     advection!(fᵗ, p, meshv, e, nx, dt)
     transpose!(f, fᵗ)
     advection!(f, p, meshx, v, nv, 0.5*dt)
     push!(nrj, 0.5*log(sum(e.*e)*dx))
  end

  nrj

end

using Plots, LinearAlgebra
using Splittings
pyplot()

nt = 600
tf = 60.0
t  = range(0.0, stop=tf, length=nt)
@time nrj = landau(tf, nt)
plot( t, nrj)
plot!(t, -0.1533*t.-5.50)