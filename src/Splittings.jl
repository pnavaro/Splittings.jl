module Splittings

export Mesh, RectMesh2D, bspline, compute_rho, compute_e

struct Mesh
    nx   :: Int
    ny   :: Int
    xmin :: Float64
    xmax :: Float64
    ymin :: Float64
    ymax :: Float64
end

"""

   2D rectangular cartesian mesh parameters

"""
struct RectMesh2D

    xmin :: Float64
    xmax :: Float64
    nx   :: Int
    dx   :: Float64
    ymin :: Float64
    ymax :: Float64
    ny   :: Int
    dy   :: Float64

    function RectMesh2D(xmin, xmax, nx, ymin, ymax, ny)
       dx = (xmax - xmin) / nx
       dy = (ymax - ymin) / ny
       new( xmin, xmax, nx, dx, ymin, ymax, ny, dy)
    end

end

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

struct UniformMesh
    left::Float64
    right::Float64
    ncells::Int32
end
delta(mesh::UniformMesh) = (mesh.right - mesh.left) / mesh.ncells
getpoints(mesh::UniformMesh) = range(mesh.left, stop=mesh.right,length=mesh.ncells+1)[1:end-1]


"""
Advection in υ
∂ f / ∂ t − E(x) ∂ f / ∂ υ  = 0
"""
function advection_v!( fᵀ, meshx::UniformMesh, meshv::UniformMesh, E, dt)

    n = meshv.ncells
    L = meshv.right - meshv.left
    k = 2π/L*[0:n÷2-1;-n÷2:-1]
    ek = exp.(-1im * dt * k * transpose(E))

    fft!(fᵀ, 1)
    fᵀ .= fᵀ .* ek
    ifft!(fᵀ, 1)

end

function advection_x!( f, meshx, meshv, dt)

    L = meshx.right - meshx.left
    m = div(meshx.ncells,2)
    k = 2*π/L * [0:1:m-1;-m:1:-1]
    k̃ = 2*π/L * [1;1:1:m-1;-m:1:-1]
    v = getpoints(meshv)
    ev = exp.(-1im*dt * k * transpose(v))

    fft!(f,1)
    f .= f .* ev
    Ek  = -1im * delta(meshv) * sum(f,dims=2) ./ k̃
    Ek[1] = 0.0
    ifft!(f,1)
    real(ifft(Ek))

end


"""
Compute charge density
ρ(x,t) = ∫ f(x,v,t) dv
"""
function compute_rho(meshv, f)
   dv = meshv.dx
   rho = dv * sum(f, dims=2)
   rho .- mean(rho)
end

"""
 compute Ex using that -ik*Ex = rho
"""
function compute_e(meshx, rho)
   nx = meshx.nx
   k =  2* pi / (meshx.xmax - meshx.xmin)
   modes = zeros(Float64, nx)
   modes .= k * vcat(0:div(nx,2)-1,-div(nx,2):-1)
   modes[1] = 1.0
   rhok = fft(rho)./modes
   real(ifft(-1im*rhok))
end




end # module
