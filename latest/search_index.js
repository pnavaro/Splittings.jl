var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Splittings.@Lie",
    "page": "Introduction",
    "title": "Splittings.@Lie",
    "category": "macro",
    "text": "@Lie( push_t, push_v )\n\nApply the first order Lie splitting\n\npush_t and push_v are two function calls with\n`dt` as argument.\n\n\n\n\n\n"
},

{
    "location": "index.html#Splittings.@Strang",
    "page": "Introduction",
    "title": "Splittings.@Strang",
    "category": "macro",
    "text": "@Strang( push_t, push_v )\n\nApply the second order Strang splitting\n\npush_t and push_v are two function calls with\n`dt` as argument.\n\n\n\n\n\n"
},

{
    "location": "index.html#Splittings.@TripleJump",
    "page": "Introduction",
    "title": "Splittings.@TripleJump",
    "category": "macro",
    "text": "@TripleJump( push_t, push_v )\n\nApply the fourth order Triple Jump splitting\n\npush_t and push_v are two function calls with\n`dt` as argument.\n\n\n\n\n\n"
},

{
    "location": "index.html#Splittings.@Order6",
    "page": "Introduction",
    "title": "Splittings.@Order6",
    "category": "macro",
    "text": "@Order6( push_t, push_v )\n\nApply the sixth order splitting\n\npush_t and push_v are two function calls with\n`dt` as argument.\n\n\n\n\n\n"
},

{
    "location": "index.html#Splittings.jl-Documentation-1",
    "page": "Introduction",
    "title": "Splittings.jl Documentation",
    "category": "section",
    "text": "Operators splitting package to solve equations of the form  fracdUdt = (T+V)Uwhere T and V are differential operators.The solution on one time step can be written U(Δt) = mathcalS_T+V U(0) The composition algorithm consists in successive solutions of the split equations fracdUdt = T U and fracdUdt = V U Alternating the two reduced solution operators mathcalS_T and mathcalS_V with adequately chosen time increments yields arbitrary  order in time for the full solution.The application of an operator splitting method to a concrete problem is done by using Julia macros:@Lie\n@Strang\n@TripleJump\n@Order6Examples of applications are provided for~:The linear pendulum problem.\nThe Vlasov equation with constant coefficients advection field.\nThe non linear Vlasov-Poisson equations in cartesian coordinates.This code is derived from Fortran and Python codes written by~:Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences).\nMichel Mehrenberger  (Aix-Marseille Université).\nEric Sonnendrucker   (Max-Planck-Institut für Plasmaphysik)."
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
    "location": "advections.html#",
    "page": "Advection functions",
    "title": "Advection functions",
    "category": "page",
    "text": ""
},

{
    "location": "advections.html#Advection-functions-1",
    "page": "Advection functions",
    "title": "Advection functions",
    "category": "section",
    "text": "CurrentModule = SplittingsModules = [Splittings]\nOrder   = [:advection!]"
},

{
    "location": "examples/vlasov-ampere.html#",
    "page": "Vlasov-Ampere",
    "title": "Vlasov-Ampere",
    "category": "page",
    "text": "EditURL = \"https://github.com/pnavaro/Splittings.jl/blob/master/examples/vlasov-ampere.jl\""
},

{
    "location": "examples/vlasov-ampere.html#Vlasov-Ampere-1",
    "page": "Vlasov-Ampere",
    "title": "Vlasov-Ampere",
    "category": "section",
    "text": "notebook,Compute Landau damping by solving Vlasov-Ampere system. fracft + υ fracfx - E(tx) fracfυ = 0fracEt = - J =  fυ dυ"
},

{
    "location": "examples/vlasov-ampere.html#Splittings.advection_x!",
    "page": "Vlasov-Ampere",
    "title": "Splittings.advection_x!",
    "category": "function",
    "text": "advection_x!( f, meshx, meshv, dt) \n\nAdvection in x and compute electric field\n\n∂ f / ∂ t − υ ∂ f / ∂ x  = 0\n\n∂E / ∂t = −J = ∫ fυ dυ\n\n\n\n\n\n"
},

{
    "location": "examples/vlasov-ampere.html#Splittings.advection_v!",
    "page": "Vlasov-Ampere",
    "title": "Splittings.advection_v!",
    "category": "function",
    "text": "advection_v!( fᵀ, meshx, meshv, E, dt) \n\nAdvection in υ\n\n∂ f / ∂ t − E(x) ∂ f / ∂ υ  = 0\n\nApply this advection on transposed `f`\n\n\n\n\n\n"
},

{
    "location": "examples/vlasov-ampere.html#Algorithm-1",
    "page": "Vlasov-Ampere",
    "title": "Algorithm",
    "category": "section",
    "text": "For each j compute discrete Fourier transform in x of (x_iυ_j) yielding f_k^n(υ_j),\nFor k  0\nCompute f^n+1_k(υ_j) = e^2iπ k υ ΔtL f_n^k(υ_j)\nCompute ρ_k^n+1 = Δ υ _j f^n+1_k(υ_j)\nCompute E^n+1_k = ρ^n+1_k L(2iπkϵ_0)\nFor k = 0 do nothing:f_n+1(υ_j) = f^n_k(υ_j) E^n+1_k = E^n_kPerform inverse discrete Fourier transform of E^n+1_k and for each j of f^n+1_k (υ_j).import Splittings: advection_x!, advection_v!, UniformMesh\nimport Splittings: @Strang, compute_rho, compute_e\nusing Plots, LinearAlgebra\npyplot()Splittings.advection_x!\nSplittings.advection_v!function push_t!( f, fᵀ, meshx, meshv, e,  dt)\n    transpose!(f,fᵀ)\n    advection_x!( f, meshx, meshv, e,  dt)\n    transpose!(fᵀ,f)\nendfunction push_v!(fᵀ, meshx, meshv, e, dt)\n    advection_v!(fᵀ, meshx, meshv, e, dt)\nendfunction vm1d( nx, nv, xmin, xmax, vmin, vmax , tf, nt)\n\n    meshx = UniformMesh(xmin, xmax, nx, endpoint=false)\n    meshv = UniformMesh(vmin, vmax, nv, endpoint=false)\n\n    x = meshx.points\n    v = meshv.points\n    ϵ, kx = 0.001, 0.5\n\n    f = zeros(Complex{Float64},(nx,nv))\n    fᵀ= zeros(Complex{Float64},(nv,nx))\n\n    f .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))\n    transpose!(fᵀ,f)\n\n    e = zeros(Complex{Float64},nx)\n\n    ρ  = compute_rho(meshv, f)\n    e .= compute_e(meshx, ρ)\n\n    nrj = Float64[]\n\n    dt = tf / nt\n\n    for i in 1:nt\n\n	push!(nrj, 0.5*log(sum(real(e).^2)*meshx.step))\n\n        @Strang(  push_v!(fᵀ, meshx, meshv, e,  dt),\n                  push_t!(f, fᵀ, meshx, meshv, e,  dt)\n               )\n\n    end\n    nrj\nendnx, nv = 32, 64\nxmin, xmax =  0., 4π\nvmin, vmax = -6., 6.\ntf = 80\nnt = 600\n\nt = range(0,stop=tf,length=nt)\nplot(t, vm1d(nx, nv, xmin, xmax, vmin, vmax, tf, nt) )\nplot!(t, -0.1533*t.-5.50)\nsavefig(\"va-plot.png\"); nothing(Image: )This page was generated using Literate.jl."
},

{
    "location": "examples/vlasov-poisson.html#",
    "page": "Vlasov-Poisson",
    "title": "Vlasov-Poisson",
    "category": "page",
    "text": "EditURL = \"https://github.com/pnavaro/Splittings.jl/blob/master/examples/vlasov-poisson.jl\""
},

{
    "location": "examples/vlasov-poisson.html#Vlasov-Poisson-1",
    "page": "Vlasov-Poisson",
    "title": "Vlasov-Poisson",
    "category": "section",
    "text": "notebook,We consider the dimensionless Vlasov-Poisson equation for one species with a neutralizing background. fracft+ v_x f + E(tx)  _v f = 0 \n - Δϕ = 1 - ρ E = -  ϕ \n ρ(tx)  =   f(txv)dvVlasov Equation - Wikipedia\nLandau damping - Wikipediausing Plots, LinearAlgebra\npyplot()\nimport Splittings:UniformMesh, BSpline\nimport Splittings:@Strang\nimport Splittingsfunction push_t!(f, meshx, v, nv, dt)\n    Splittings.advection!(f, meshx, v, nv, dt, BSpline(5))\nendfunction push_v!(f, fᵗ, meshx, meshv, nrj, dt)\n    rho = Splittings.compute_rho(meshv, f)\n    e   = Splittings.compute_e(meshx, rho)\n    push!(nrj, 0.5*log(sum(e.*e)*meshx.step))\n    transpose!(fᵗ, f)\n    Splittings.advection!(fᵗ, meshv, e, meshx.length, dt, BSpline(5))\n    transpose!(f, fᵗ)\nendfunction landau(tf, nt)\n\n  nx, nv = 32, 64\n  xmin, xmax = 0.0, 4π\n  vmin, vmax = -6., 6.\n  meshx = UniformMesh(xmin, xmax, nx; endpoint=false)\n  meshv = UniformMesh(vmin, vmax, nv; endpoint=false)\n  x = meshx.points\n  v = meshv.points\n  dx = meshx.step\n\n  ϵ, kx = 0.001, 0.5\n  f = zeros(Complex{Float64},(nx,nv))\n  f .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) * transpose(exp.(-0.5*v.^2))\n  fᵗ = zeros(Complex{Float64},(nv,nx))\n\n  dt = tf / nt\n\n  nrj = Float64[]\n\n  for it in 1:nt\n      @Strang( push_t!(f, meshx, v, nv, dt),\n               push_v!(f, fᵗ, meshx, meshv, nrj, dt))\n  end\n\n  nrj\n\nendnt = 600\ntf = 60.0\nt  = range(0.0, stop=tf, length=nt)\n@time nrj = landau(tf, nt)\nplot( t, nrj)\nplot!(t, -0.1533*t.-5.50)\nsavefig(\"landau-plot.png\"); nothing # hide(Image: )This page was generated using Literate.jl."
},

{
    "location": "examples/bump_on_tail.html#",
    "page": "Bump On Tail",
    "title": "Bump On Tail",
    "category": "page",
    "text": "EditURL = \"https://github.com/pnavaro/Splittings.jl/blob/master/examples/bump_on_tail.jl\""
},

{
    "location": "examples/bump_on_tail.html#Bump-On-Tail-1",
    "page": "Bump On Tail",
    "title": "Bump On Tail",
    "category": "section",
    "text": "notebook,import Splittings: advection!, UniformMesh\nimport Splittings: compute_rho, compute_e\nusing Plots\nusing LaTeXStrings\n\npyplot()function vlasov_poisson(meshx  :: UniformMesh,\n                        meshv  :: UniformMesh,\n                        f      :: Array{Float64,2},\n                        nstep  :: Int64,\n                        dt     :: Float64)\n\n    x = meshx.points\n    v = meshv.points\n\n    nrj = Float64[]\n    for istep in 1:nstep\n        advection!( f, meshx, v, 0.5dt, axis=1)\n        rho = compute_rho(meshv, f)\n        e   = compute_e(meshx, rho)\n        advection!( f, meshv, e, dt, axis=2)\n        advection!( f, meshx, v, 0.5dt, axis=1)\n        push!(nrj, 0.5*log(sum(e.*e)*meshx.step))\n    end\n    nrj\n\nendα = 0.03\nkx  = 0.3\nxmin, xmax = 0.0, 2π / kx\nnx, nv = 32, 64\nvmin, vmax = -9., 9.\nmeshx = UniformMesh(xmin, xmax, nx)\nmeshv = UniformMesh(vmin, vmax, nv)\nf = zeros(Float64,(meshx.length,meshv.length))\nfor (i,x) in enumerate(meshx.points), (j,v) in enumerate(meshv.points)\n    f[i,j]  = (1.0+α*cos(kx*x)) / (10*sqrt(2π)) * (9*exp(-0.5*v^2)+2*exp(-2*(v-4.5)^2))\nendnstep = 500\nt = range(0.0, stop=50.0, length=nstep)\ndt = t[2]\n@elapsed nrj = vlasov_poisson( meshx, meshv, f, nstep, dt)plot(t, nrj, label=L\"\\frac{1}{2} \\log(∫e²dx)\")\nsavefig(\"bot-plot.png\"); nothing # hide(Image: )This page was generated using Literate.jl."
},

{
    "location": "examples/rotation2d_bsl.html#",
    "page": "Rotation 2D",
    "title": "Rotation 2D",
    "category": "page",
    "text": "EditURL = \"https://github.com/pnavaro/Splittings.jl/blob/master/examples/rotation2d_bsl.jl\""
},

{
    "location": "examples/rotation2d_bsl.html#Rotation-of-a-gaussian-distribution-1",
    "page": "Rotation 2D",
    "title": "Rotation of a gaussian distribution",
    "category": "section",
    "text": "#md # notebook,    fracdfdt +  (y fracdfdx - x fracdfdy) = 0import Splittings: advection!, UniformMesh\nimport Splittings: @Strang\nusing Plots\npyplot()function with_bsl(tf::Float64, nt::Int)\n\n   nx, ny = 32, 64\n   meshx = UniformMesh(-π, π, nx)\n   meshy = UniformMesh(-π, π, ny)\n   x = meshx.points\n   y = meshy.points\n\n   dt = tf/nt\n\n   f = zeros(Float64,(nx,ny))\n\n   for (i, xp) in enumerate(x), (j, yp) in enumerate(y)\n       xn = cos(tf)*xp - sin(tf)*yp\n       yn = sin(tf)*xp + cos(tf)*yp\n       f[i,j] = exp(-(xn-1)*(xn-1)/0.2)*exp(-(yn-1)*(yn-1)/0.2)\n   end\n\n   anim = @animate for n=1:nt\n\n      @Strang(advection!( f,  meshx,  y, tan(dt), axis=1),\n              advection!( f,  meshy, -x, sin(dt), axis=2))\n\n      surface(f)\n\n   end\n\n   gif(anim, \"rotanim.gif\", fps=15); nothing #hide\n\nend@time f = with_bsl( 2π, 20)(Image: )This page was generated using Literate.jl."
},

{
    "location": "examples/vlasov-hmf.html#",
    "page": "Vlasov-HMF",
    "title": "Vlasov-HMF",
    "category": "page",
    "text": "EditURL = \"https://github.com/pnavaro/Splittings.jl/blob/master/examples/vlasov-hmf.jl\""
},

{
    "location": "examples/vlasov-hmf.html#Vlasov-HMF-1",
    "page": "Vlasov-HMF",
    "title": "Vlasov-HMF",
    "category": "section",
    "text": "notebook,using LinearAlgebra, QuadGK, Roots, FFTW\nusing Splittings\nusing Plots\npyplot()\" Compute M₀ by solving F(m) = 0 \"\nfunction mag(β, mass)\n\n    F(m) = begin\n        g(x, n, m) = (1 / π) * (exp(β * m * cos(x)) * cos(n * x))\n        bessel0(x) = g(x, 0, m)\n        bessel1(x) = g(x, 1, m)\n        mass * quadgk(bessel1, 0, π)[1] / quadgk(bessel0, 0, π)[1] - m\n    end\n\n    find_zero(F, (0, mass))\nendfunction Norm(f::Array{Float64,2}, dx, dv)\n   return dx * sum(dv * sum(real(f), dims=1))\nend\"\"\"\n    Compute the electric hamiltonian mean field from a\n    2D distribution function\n\"\"\"\nfunction hmf_poisson!(fᵗ::Array{Complex{Float64},2},\n        meshx::UniformMesh,\n        meshv::UniformMesh,\n        ex::Array{Float64})\n\n    nx = meshx.length\n    rho = meshv.step .* vec(sum(fᵗ, dims=1))\n    kernel = zeros(Float64, nx)\n    k = π / (meshx.stop - meshx.start)\n    kernel[2] = k\n    ex .= real(ifft(1im * fft(rho) .* kernel * 4π ))\n\nendfunction bsl_advection!(f::Array{Complex{Float64},2},\n                        meshx::UniformMesh,\n                        meshv::UniformMesh,\n                        v::Array{Float64,1},\n                        dt;\n                        spline_degree=3)\n\n    fft!(f,1)\n    @simd for j in 1:meshv.length\n        alpha = v[j] * dt\n        @inbounds f[:,j] .= Splittings.interpolate(spline_degree, f[:,j],\n            meshx.step, alpha)\n    end\n    ifft!(f,1)\nendfunction push_v!(f, fᵗ, meshx, meshv, ex, dt)\n    transpose!(fᵗ, f)\n    hmf_poisson!(fᵗ, meshx, meshv, ex)\n    bsl_advection!(fᵗ, meshv, meshx, ex, dt)\n    transpose!(f, fᵗ)\nendfunction vlasov_hmf_gauss(nbiter = 10000, dt = 0.1)\n\n    mass = 1.0\n    T = 0.1\n    meshx = UniformMesh(-π, π, 64)\n    meshv = UniformMesh(-8, 8, 64)\n\n    nx, dx = meshx.length, meshx.step\n    nv, dv = meshv.length, meshv.step\n    x, v = meshx.points, meshv.points\n    X = repeat(x,1,nv)\n    V = repeat(v,1,nx)\'\n    ϵ = 0.1\n\n    b = 1 / T\n    m = mag(b, mass)\n\n    w   = sqrt(m)\n    f   = zeros(Complex{Float64}, (nx,nv))\n    fᵗ  = zeros(Complex{Float64}, (nv,nx))\n\n    f  .= exp.(-b .* ((V.^2 / 2) - m * cos.(X)))\n    a   = mass / Norm(real(f), dx, dv)\n    @.  f =  a * exp(-b * (((V^2) / 2) - m * cos(X))) * (1 + ϵ * cos(X))\n\n    ex = zeros(Float64,nx)\n    hmf_poisson!(f, meshx, meshv, ex )\n    test = copy(f)\n    T = Float64[]\n    for n in 1:nbiter\n\n        gamma1 = Norm(real(f) .* cos.(X), dx, dv)\n        push!(T,gamma1)\n\n        @Strang(\n            bsl_advection!(f, meshx, meshv, v, dt),\n            push_v!(f, fᵗ, meshx, meshv, ex, dt)\n        )\n\n    end\n\n    #Substracting from gamma its long time average\n\n    Gamma1 = Norm(real(f) .*cos.(X), dx, dv)\n    T .= T .- Gamma1\n\n    range(0., stop=nbiter*deltat, length=nbiter), abs.(T)\n\nendnbiter = 2000\ndeltat = 0.1\n@time t, T = vlasov_hmf_gauss(nbiter, deltat);\nplot(t, log.(T), xlabel = \"t\", ylabel = \"|C[f](t)-C[f][T]|\")\nsavefig(\"vlasov-hmf-plot.png\"); nothing # hide(Image: png)This page was generated using Literate.jl."
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

{
    "location": "contents.html#",
    "page": "Contents",
    "title": "Contents",
    "category": "page",
    "text": ""
},

{
    "location": "contents.html#Contents-1",
    "page": "Contents",
    "title": "Contents",
    "category": "section",
    "text": ""
},

{
    "location": "contents.html#Index-1",
    "page": "Contents",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
