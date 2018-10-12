var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Splittings.jl-Documentation-1",
    "page": "Introduction",
    "title": "Splittings.jl Documentation",
    "category": "section",
    "text": "Operators splitting package to solve equations of the form dUdt = (T+V)U,where T and  V are two differential operators by solving successively the simpler equationsThe application of an operator splitting method to a concrete problem is done by using Julia macros."
},

{
    "location": "index.html#Examples-of-applications-are-provided-for-1",
    "page": "Introduction",
    "title": "Examples of applications are provided for",
    "category": "section",
    "text": "The linear pendulum problem.\nThe Vlasov equation with constant coefficients advection field.\nThe non linear Vlasov-Poisson equations in cartesian coordinates.ReferencesE. Hairer, C. Lubich, G. Wanner, Geometrical numerical integration, Springer 2006This code is derived from Fortran and Python codes written by - Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences).\n- Michel Mehrenberger  (Aix-Marseille Université).\n- Eric Sonnendrucker   (Max-Planck-Institut für Plasmaphysik).CurrentModule = SplittingsSemi-Lagrangian method"
},

{
    "location": "index.html#Examples-1",
    "page": "Introduction",
    "title": "Examples",
    "category": "section",
    "text": "Vlasov-Poisson\nVlasov-Ampere\nBump On Tail"
},

{
    "location": "index.html#Splittings.compute_rho-Tuple{Any,Any}",
    "page": "Introduction",
    "title": "Splittings.compute_rho",
    "category": "method",
    "text": "compute_rho( mesh, f)\n\nCompute charge density\n\nρ(x,t) = ∫ f(x,v,t) dv\n\nreturn ρ - ρ̄ if neutralized=true\n\n\n\n\n\n"
},

{
    "location": "index.html#Splittings.compute_e-Tuple{Any,Any}",
    "page": "Introduction",
    "title": "Splittings.compute_e",
    "category": "method",
    "text": "compute_e( mesh, ρ)\n\nCompute 1d electric field using that -ik * e = ρ\n\n\n\n\n\n"
},

{
    "location": "index.html#Functions-1",
    "page": "Introduction",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [Splittings]\nOrder   = [:advection!, :UniformMesh]compute_rho( meshv, f)compute_e( meshx, rho)"
},

{
    "location": "index.html#Index-1",
    "page": "Introduction",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "vlasov-ampere.html#",
    "page": "Vlasov-Ampere",
    "title": "Vlasov-Ampere",
    "category": "page",
    "text": ""
},

{
    "location": "vlasov-ampere.html#Vlasov-Ampere-1",
    "page": "Vlasov-Ampere",
    "title": "Vlasov-Ampere",
    "category": "section",
    "text": "Compute Landau damping by solving Vlasov-Ampere system. fracft + υ fracfx\n - E(tx) fracfυ = 0\n \nmath \\frac{∂E}{∂t} = - J = ∫ fυ dυ\n## Algorithm\n\n - For each ``j`` compute discrete Fourier transform in ``x`` of\n   ``(x_i,υ_j)`` yielding ``f_k^n(υ_j)``,\n - For `` k ≂̸ 0 ``\n\n     - Compute\nmath f^{n+1}k(υj) = e^{−2iπ k υ Δt/L} fn^k(υj), \n     - Compute\nmath ρk^{n+1} = Δ υ ∑j􏰄 f^{n+1}k(υj), \n     - Compute\nmath E^{n+1}k = ρ^{n+1}k L/(2iπkϵ_0), \n - For ``k = 0`` do nothing:\nmath f{n+1}(υj) = f^nk(υj), E^{n+1}k = E^nk. \n - Perform inverse discrete Fourier transform of ``E^{n+1}_k`` and for each\n   ``j`` of ``f^{n+1}_k (υ_j)``.\n\n## Variables\n\n- `nx::Integer`: the number of points along x axis.\n- `nv::Integer`: the number of points along ⁠υ axis.\n- `xmin::Float`: the origin of mesh along x axis.\n- `xmax::Float`: the end of mesh along x axis.\n- `vmin::Float`: the origin of mesh along υ axis.\n- `vmax::Float`: the end of mesh along υ axis.\n- `tf::Float`  : the final time of the simulation.\n- `nt::Integer`: the number of time steps.\n@example import Splittings: advectionx!, advectionv!, UniformMesh import Splittings: @Strang, computerho, computee using Plots, LinearAlgebra pyplot()function pusht!( f, fᵀ, meshx, meshv, e,  dt)     transpose!(f,fᵀ)     advectionx!( f, meshx, meshv, e,  dt)     transpose!(fᵀ,f) endfunction pushv!(fᵀ, meshx, meshv, e, dt)     advectionv!(fᵀ, meshx, meshv, e, dt) endfunction vm1d( nx, nv, xmin, xmax, vmin, vmax , tf, nt)meshx = UniformMesh(xmin, xmax, nx, endpoint=false)\nmeshv = UniformMesh(vmin, vmax, nv, endpoint=false)\n\nx = meshx.x\nv = meshv.x\nϵ, kx = 0.001, 0.5\n\nf = zeros(Complex{Float64},(nx,nv))\nfᵀ= zeros(Complex{Float64},(nv,nx))\n\nf .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))\ntranspose!(fᵀ,f)\n\ne = zeros(Complex{Float64},nx)\n\nρ  = compute_rho(meshv, f)\ne .= compute_e(meshx, ρ)\n\nnrj = Float64[]\n\ndt = tf / nt\n\nfor i in 1:nt\n\npush!(nrj, 0.5*log(sum(real(e).^2)*meshx.dx))\n\n    @Strang(  push_v!(fᵀ, meshx, meshv, e,  dt),\n              push_t!(f, fᵀ, meshx, meshv, e,  dt)\n           )\n              \nend\nnrjendnx, nv = 64, 128 xmin, xmax =  0., 4π vmin, vmax = -6., 6. tf = 80 nt = 600t = range(0,stop=tf,length=nt) plot(t, vm1d(nx, nv, xmin, xmax, vmin, vmax, tf, nt) ) plot!(t, -0.1533*t.-5.50) savefig(\"va-plot.png\"); nothing # hide ```(Image: )"
},

{
    "location": "vlasov-poisson.html#",
    "page": "Vlasov-Poisson",
    "title": "Vlasov-Poisson",
    "category": "page",
    "text": ""
},

{
    "location": "vlasov-poisson.html#Vlasov-Poisson-1",
    "page": "Vlasov-Poisson",
    "title": "Vlasov-Poisson",
    "category": "section",
    "text": "\nusing Plots, LinearAlgebra\nusing Splittings\n\n\n\"\"\"\n    landau(tf, nt)\n\n Compute Landau damping by solving Vlasov-Poisson system in 1D1V space.\n \n ## Semi-Lagrangian method\n\n Let us consider an abstract scalar advection equation of the form\n \n ``\\\\frac{∂f}{∂t}+ a(x, t) ⋅ ∇f = 0.``\n \n The characteristic curves associated to this equation are the solutions of \n the ordinary differential equations\n \n ``\\\\frac{dX}{dt} = a(X(t), t)``\n\n We shall denote by ``X(t, x, s)`` the unique solution of this equation \n associated to the initial condition ``X(s) = x``.\n\n The classical semi-Lagrangian method is based on a backtracking of \n characteristics. Two steps are needed to update the distribution function \n ``f^{n+1}`` at ``t^{n+1}`` from its value ``f^n`` at time ``t^n`` :\n\n 1. For each grid point ``x_i`` compute ``X(t^n; x_i, t^{n+1})`` the value \n    of the characteristic at ``t^n`` which takes the value ``x_i`` at \n    ``t^{n+1}``.\n 2. As the distribution solution of first equation verifies\n    ``f^{n+1}(x_i) = f^n(X(t^n; x_i, t^{n+1})),``\n    we obtain the desired value of ``f^{n+1}(x_i)`` by computing \n    ``f^n(X(t^n;x_i,t^{n+1})`` by interpolation as ``X(t^n; x_i, t^{n+1})`` \n    is in general not a grid point.\n\n *[Eric Sonnendrücker - Numerical methods for the Vlasov equations](http://www-m16.ma.tum.de/foswiki/pub/M16/Allgemeines/NumMethVlasov/Num-Meth-Vlasov-Notes.pdf)*\n\n\n Vlasov-Poisson equation\n -----------------------\n\n We consider the dimensionless Vlasov-Poisson equation for one species\n with a neutralizing background.\n\n ``\n \\\\frac{∂f}{∂t}+ v⋅∇_x f + E(t,x) ⋅ ∇_v f = 0, \\\\\n - Δϕ = 1 - ρ, E = - ∇ ϕ \\\\\n ρ(t,x)  =  ∫ f(t,x,v)dv.\n ``\n\n - [Vlasov Equation - Wikipedia](https://en.wikipedia.org/wiki/Vlasov_equation)\n - [Landau damping - Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)\n\n# Examplesjulia nt = 600 tf = 60.0 t  = range(0.0, stop=tf, length=nt) @time nrj = landau(tf, nt) plot( t, nrj) plot!(t, -0.1533*t.-5.50)\n\"\"\"\nfunction landau(tf, nt)\n\n  # Set grid\n  p = 3\n  nx, nv = 64, 128\n  xmin, xmax = 0.0, 4π\n  vmin, vmax = -6., 6.\n  meshx = UniformMesh(xmin, xmax, nx; endpoint=false)\n  meshv = UniformMesh(vmin, vmax, nv; endpoint=false)\n  x = meshx.x\n  v = meshv.x\n  dx = meshx.dx\n\n  # Initialize distribution function for Landau damping.\n\n  eps, kx = 0.001, 0.5\n  f = zeros(Complex{Float64},(nx,nv))\n  f .= (1.0.+eps*cos.(kx*x))/sqrt(2.0*pi) * transpose(exp.(-0.5*v.^2))\n  fᵗ = zeros(Complex{Float64},(nv,nx))\n\n  # Set time domain\n  dt = tf / nt\n\n  # Run simulation\n  nrj = Float64[]\n\n  for it in 1:nt\n     Splittings.advection!(f, p, meshx, v, nv, 0.5*dt)\n     rho = Splittings.compute_rho(meshv, f)\n     e   = Splittings.compute_e(meshx, rho)\n     transpose!(fᵗ, f)\n     Splittings.advection!(fᵗ, p, meshv, e, nx, dt)\n     transpose!(f, fᵗ)\n     Splittings.advection!(f, p, meshx, v, nv, 0.5*dt)\n     push!(nrj, 0.5*log(sum(e.*e)*dx))\n  end\n\n  nrj\n\nend\n\nusing Plots, LinearAlgebra\npyplot()\n\nnt = 600\ntf = 60.0\nt  = range(0.0, stop=tf, length=nt)\n@time nrj = landau(tf, nt)\nplot( t, nrj)\nplot!(t, -0.1533*t.-5.50)\nsavefig(\"landau-plot.png\"); nothing # hide(Image: )"
},

{
    "location": "bump_on_tail.html#",
    "page": "Bump On Tail",
    "title": "Bump On Tail",
    "category": "page",
    "text": ""
},

{
    "location": "bump_on_tail.html#Bump-On-Tail-1",
    "page": "Bump On Tail",
    "title": "Bump On Tail",
    "category": "section",
    "text": "\nusing Splittings\nusing Plots\n\npyplot()\n\n\nfunction vlasov_poisson(mesh  :: RectMesh1D1V, \n                        f     :: Array{Float64,2}, \n                        nstep :: Int64, \n                        dt    :: Float64)\n    \n    nrj = Float64[]\n    advection!( f, mesh, 0.5dt, axis=1)\n    for istep in 1:nstep\n        rho = compute_rho(mesh, f)\n        e   = compute_e(mesh, rho)\n        advection!( f, mesh, e, dt)\n        advection!( f, mesh, dt)\n        push!(nrj, 0.5*log(sum(e.*e)*mesh.dx))\n    end        \n    nrj\n    \nend\n\nα = 0.03\nkx  = 0.3\nxmin, xmax = 0.0, 2π / kx\nnx, nv = 512, 512\nvmin, vmax = -9., 9.\nmesh = RectMesh1D1V(xmin, xmax, nx, vmin, vmax, nv)\nf = zeros(Float64,(mesh.nx,mesh.nv))           \nfor (i,x) in enumerate(mesh.x), (j,v) in enumerate(mesh.v)\n     f[i,j]  = (1.0+α*cos(kx*x)) / (10*sqrt(2π)) * (9*exp(-0.5*v^2)+2*exp(-2*(v-4.5)^2))\nend\n\n\nnstep = 500\nt = range(0.0, stop=50.0, length=nstep)\ndt = t[2]\n@elapsed nrj = vlasov_poisson( mesh, f, nstep, dt)\n\nplot(t, nrj, label=\"``\\\\frac{1}{2}log(∫e²dx)``\")\nsavefig(\"bot-plot.png\"); nothing # hide(Image: )"
},

{
    "location": "rotation2d_bsl.html#",
    "page": "Rotation 2D",
    "title": "Rotation 2D",
    "category": "page",
    "text": ""
},

{
    "location": "rotation2d_bsl.html#Rotation-of-a-gaussian-distribution-1",
    "page": "Rotation 2D",
    "title": "Rotation of a gaussian distribution",
    "category": "section",
    "text": "    fracdfdt +  (y fracdfdx - x fracdfdy) = 0import Splittings: advection!, UniformMesh\nimport Splittings: @Strang\nusing Plots\npyplot(leg=false, ticks=nothing)\n\nfunction with_bsl(tf::Float64, nt::Int)\n\n   nx, ny = 64, 64\n   meshx = UniformMesh(-π, π, nx)\n   meshy = UniformMesh(-π, π, ny)\n   x = meshx.x\n   y = meshy.x\n\n   dt = tf/nt\n\n   f = zeros(Float64,(nx,ny))\n\n   for (i, xp) in enumerate(x), (j, yp) in enumerate(y)\n       xn = cos(tf)*xp - sin(tf)*yp\n       yn = sin(tf)*xp + cos(tf)*yp\n       f[i,j] = exp(-(xn-1)*(xn-1)/0.2)*exp(-(yn-1)*(yn-1)/0.2)\n   end\n\n   anim = @animate for n=1:nt\n       \n      @Strang(advection!( f,  meshx,  y, tan(dt), axis=1),\n              advection!( f,  meshy, -x, sin(dt), axis=2))\n\n      surface(f)\n                                      \n   end\n\n   gif(anim, \"rotanim.gif\", fps=15)\n\nend\n\nf = with_bsl( 10π, 100)(Image: )"
},

{
    "location": "bsl.html#",
    "page": "Semi-Lagrangian",
    "title": "Semi-Lagrangian",
    "category": "page",
    "text": ""
},

{
    "location": "bsl.html#Semi-Lagrangian-method-1",
    "page": "Semi-Lagrangian",
    "title": "Semi-Lagrangian method",
    "category": "section",
    "text": "Let us consider an abstract scalar advection equation of the form$ \\frac{∂f}{∂t}+ a(x, t) ⋅ ∇f = 0. $The characteristic curves associated to this equation are the solutions of  the ordinary differential equations$ \\frac{dX}{dt} = a(X(t), t) $We shall denote by X(t x s) the unique solution of this equation  associated to the initial condition X(s) = x.The classical semi-Lagrangian method is based on a backtracking of  characteristics. Two steps are needed to update the distribution function  f^n+1 at t^n+1 from its value f^n at time t^n :For each grid point x_i compute X(t^n x_i t^n+1) the value  of the characteristic at t^n which takes the value x_i at  t^n+1.\nAs the distribution solution of first equation verifies f^n+1(x_i) = f^n(X(t^n x_i t^n+1)) we obtain the desired value of f^n+1(x_i) by computing  f^n(X(t^nx_it^n+1) by interpolation as X(t^n x_i t^n+1)  is in general not a grid point.Eric Sonnendrücker - Numerical methods for the Vlasov equationsCurrentModule = Splittings"
},

{
    "location": "bsl.html#Functions-1",
    "page": "Semi-Lagrangian",
    "title": "Functions",
    "category": "section",
    "text": "advection!(f, p, mesh, v, nv, dt)"
},

{
    "location": "bsl.html#Index-1",
    "page": "Semi-Lagrangian",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "contributing.html#",
    "page": "How to Contribute",
    "title": "How to Contribute",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#How-to-Contribute-1",
    "page": "How to Contribute",
    "title": "How to Contribute",
    "category": "section",
    "text": "Here\'s an outline of the workflow you should use if you want to make  contributions to Splittings.Fork Splittings\nMake a new branch on your fork, named after whatever changes you\'ll be making\nApply your code changes to the branch on your fork\nWhen you\'re done, submit a PR to Splittings to merge your fork into  Splittings\'s master branch."
},

]}
