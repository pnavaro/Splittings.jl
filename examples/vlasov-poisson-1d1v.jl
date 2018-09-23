#
#
# Semi-Lagrangian method
# ----------------------
#
# Let us consider an abstract scalar advection equation of the form
# $$
# \frac{\partial f}{\partial t}+ a(x, t) \cdot \nabla f = 0.
# $$
# The characteristic curves associated to this equation are the solutions of the ordinary differential equations
# $$
# \frac{dX}{dt} = a(X(t), t)
# $$
# We shall denote by $X(t, x, s)$ the unique solution of this equation associated to the initial condition $X(s) = x$.
#
# The classical semi-Lagrangian method is based on a backtracking of characteristics. Two steps are needed to update the distribution function $f^{n+1}$ at $t^{n+1}$ from its value $f^n$ at time $t^n$ :
#
# 1. For each grid point $x_i$ compute $X(t^n; x_i, t^{n+1})$ the value of the characteristic at $t^n$ which takes the value $x_i$ at $t^{n+1}$.
# 2. As the distribution solution of first equation verifies
# $$f^{n+1}(x_i) = f^n(X(t^n; x_i, t^{n+1})),$$
# we obtain the desired value of $f^{n+1}(x_i)$ by computing $f^n(X(t^n;x_i,t^{n+1})$ by interpolation as $X(t^n; x_i, t^{n+1})$ is in general not a grid point.
#
# *[Eric Sonnendrücker - Numerical methods for the Vlasov equations](http://www-m16.ma.tum.de/foswiki/pub/M16/Allgemeines/NumMethVlasov/Num-Meth-Vlasov-Notes.pdf)*
#
#



using Plots, FFTW, LinearAlgebra, Statistics
using Splittings


"""
[De Boor's Algorithm - Wikipedia](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm)

Return the value at x in [0,1[ of the B-spline with
integer nodes of degree p with support starting at j.
Implemented recursively using the de Boor's recursion formula
"""
function bspline(p::Int, j::Int, x::Float64)
   if p == 0
       if j == 0
           return 1.0
       else
           return 0.0
       end
   else
       w = (x - j) / p
       w1 = (x - j - 1) / p
   end
   return (w * bspline(p - 1, j, x) + (1 - w1) * bspline(p - 1, j + 1, x))
end

"""
1D uniform mesh data
"""
struct UniformMesh
   xmin  :: Float64
   xmax  :: Float64
   nx    :: Int
   dx    :: Float64
   x     :: Vector{Float64}
   function UniformMesh(xmin, xmax, nx)
      dx = (xmax - xmin) / nx
      x  = range(xmin, stop=xmax, length=nx+1)[1:end-1]
      new( xmin, xmax, nx, dx, x)
   end
end

function advection!(f, p, mesh, v, nv, dt)

   nx = mesh.nx
   dx = mesh.dx
   modes = [2π * i / nx for i in 0:nx-1]
   # compute eigenvalues of degree p b-spline matrix
   eig_bspl = zeros(Float64, nx)
   eig_bspl .= bspline(p, -div(p+1,2), 0.0)
   for i in 1:div(p+1,2)-1
      eig_bspl .+= bspline(p, i - div(p+1,2), 0.0) * 2 .* cos.(i * modes)
   end
   eigalpha = zeros(Complex{Float64}, nx)

   fft!(f,1)

   for j in 1:nv
      alpha = dt * v[j] / dx

      # compute eigenvalues of cubic splines evaluated at displaced points
      ishift = floor(-alpha)
      beta   = -ishift - alpha
      fill!(eigalpha,0.0im)
      for i in -div(p-1,2):div(p+1,2)
         eigalpha .+= (bspline(p, i-div(p+1,2), beta)
                        .* exp.((ishift+i) * 1im .* modes))
      end

      # compute interpolating spline using fft and properties of circulant matrices

      f[:,j] .*= eigalpha ./ eig_bspl

   end

   ifft!(f,1)

end

#
# Vlasov-Poisson equation
# -----------------------
#
# We consider the dimensionless Vlasov-Poisson equation for one species
# with a neutralizing background.
#
# $$
# \frac{\partial f}{\partial t}+ v\cdot \nabla_x f + E(t,x) \cdot \nabla_v f = 0, \\
# - \Delta \phi = 1 - \rho, E = - \nabla \phi \\
# \rho(t,x)  =  \int f(t,x,v)dv.
# $$
#
# - [Vlasov Equation - Wikipedia](https://en.wikipedia.org/wiki/Vlasov_equation)
#
#


"""
Landau Damping

[Landau damping - Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)

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

  # Create Vlasov-Poisson simulation

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

nt = 600
tf = 60.0
t  = range(0.0, stop=tf, length=nt)
@time nrj = landau(tf, nt)
plot( t, nrj)
plot!(t, -0.1533*t.-5.50)
