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
    "text": "CurrentModule = SplittingsSemi-Lagrangian method"
},

{
    "location": "index.html#Examples-1",
    "page": "Introduction",
    "title": "Examples",
    "category": "section",
    "text": "Vlasov Poisson\nLandau Damping"
},

{
    "location": "index.html#Splittings.compute_rho-Tuple{Any,Any}",
    "page": "Introduction",
    "title": "Splittings.compute_rho",
    "category": "method",
    "text": "Compute charge density    ρ(x,t) = ∫ f(x,v,t) dv\n\n\n\n\n\n"
},

{
    "location": "index.html#Splittings.compute_e-Tuple{Any,Any}",
    "page": "Introduction",
    "title": "Splittings.compute_e",
    "category": "method",
    "text": "Compute Ex using that -ik*Ex = rho\n\n\n\n\n\n"
},

{
    "location": "index.html#Functions-1",
    "page": "Introduction",
    "title": "Functions",
    "category": "section",
    "text": "compute_rho( meshv, f)compute_e( meshx, rho)"
},

{
    "location": "index.html#Index-1",
    "page": "Introduction",
    "title": "Index",
    "category": "section",
    "text": ""
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
