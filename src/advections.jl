using FFTW, LinearAlgebra

"""

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
   ( w       * bspline(p - 1, j    , x) +
    (1 - w1) * bspline(p - 1, j + 1, x))
end

"""

   Compute the interpolating spline of degree p of odd
   degree of a function f on a periodic uniform mesh, at
   all points xi-alpha

"""
function interpolate(p::Int, f::Vector{Float64}, delta::Float64, alpha::Float64)

   n = size(f)[1]
   modes = 2 * pi * (0:n-1) / n
   eig_bspl = zeros(Complex{Float64},n)
   eig_bspl .= bspline(p, -div(p+1,2), 0.0)
   for j in 1:div(p+1,2)-1
      eig_bspl .+= (bspline(p, j-div(p+1,2), 0.0)
         * 2 * cos.(j * modes))
   end
   ishift = floor(- alpha / delta)
   beta = - ishift - alpha / delta
   eigalpha = zeros(Complex{Float64},n)
   for j in -div(p-1,2):div(p+1,2)
      eigalpha .+= (bspline(p, j-div(p+1,2), beta)
         .* exp.((ishift + j) * 1im .* modes))
   end

   real(ifft(fft(f) .* eigalpha ./ eig_bspl))

end

function advection_x!(mesh::RectMesh2D, f::Array{Float64,2},
                      v::Any, dt::Float64)
    for j in 1:mesh.ny
        alpha = v[j] * dt
        f[:,j] .= interpolate(3, f[:,j], mesh.dx, alpha)
    end
end

function advection_y!(mesh::RectMesh2D, f::Array{Float64,2},
                      v::Any, dt::Float64)
    for i in 1:mesh.nx
        alpha = v[i] * dt
        f[i,:] .= interpolate(3, f[i,:], mesh.dy,  alpha)
    end
end

"""

   Advection in υ
   ∂ f / ∂ t − E(x) ∂ f / ∂ υ  = 0

"""
function advection_v!( fᵀ, meshx::UniformMesh, meshv::UniformMesh, E, dt)

    n = meshv.nx
    L = meshv.xmax - meshv.xmin
    k = 2π/L*[0:n÷2-1;-n÷2:-1]
    ek = exp.(-1im * dt * k * transpose(E))

    fft!(fᵀ, 1)
    fᵀ .= fᵀ .* ek
    ifft!(fᵀ, 1)

end

function advection_x!( f, meshx::UniformMesh, meshv::UniformMesh, dt)

    L = meshx.xmax - meshx.xmin
    m = div(meshx.nx,2)
    k = 2*π/L * [0:1:m-1;-m:1:-1]
    k̃ = 2*π/L * [1;1:1:m-1;-m:1:-1]
    v = meshv.x
    ev = exp.(-1im*dt * k * transpose(v))

    fft!(f,1)
    f .= f .* ev
    Ek  = -1im * meshv.dx * sum(f,dims=2) ./ k̃
    Ek[1] = 0.0
    ifft!(f,1)
    real(ifft(Ek))

end

function advection!(f, p, mesh, v, nv, dt)

   nx = mesh.nx
   dx = mesh.dx
   modes = [2π * i / nx for i in 0:nx-1]
   # compute eigenvalues of degree p b-spline matrix
   eig_bspl = zeros(Float64, nx)
   eig_bspl .= Splittings.bspline(p, -div(p+1,2), 0.0)
   for i in 1:div(p+1,2)-1
      eig_bspl .+= Splittings.bspline(p, i - div(p+1,2), 0.0) * 2 .* cos.(i * modes)
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
         eigalpha .+= (Splittings.bspline(p, i-div(p+1,2), beta)
                        .* exp.((ishift+i) * 1im .* modes))
      end

      # compute interpolating spline using fft and properties of circulant matrices

      f[:,j] .*= eigalpha ./ eig_bspl

   end

   ifft!(f,1)

end
