using Plots

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

   Exact solution of f after rotation during time tf

"""
function exact(tf::Float64, nt::Int, mesh::RectMesh2D)

    dt = tf/nt
    
    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    x  = linspace(xmin, xmax, nx+1)[1:end-1]

    ny = mesh.ny
    ymin, ymax = mesh.ymin, mesh.ymax
    y  = linspace(ymin, ymax, ny+1)[1:end-1]

    f = zeros(Float64,(nx,ny))
    for (i, xx) in enumerate(x), (j, yy) in enumerate(y)
        xn=cos(tf)*xx-sin(tf)*yy
        yn=sin(tf)*xx+cos(tf)*yy
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn-1)*(yn-1)/0.1)
    end

    f

end

"""
 
   Compute L1 error

"""
function error1(f, f_exact)
    maximum(abs.(f - f_exact))
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

function with_bsl(tf::Float64, nt::Int, mesh::RectMesh2D)

   dt = tf/nt
   nx = mesh.nx
   xmin, xmax = mesh.xmin, mesh.xmax
   x  = linspace(xmin, xmax, nx+1)[1:end-1]

   ny = mesh.ny
   ymin, ymax = mesh.ymin, mesh.ymax
   y  = linspace(ymin, ymax, ny+1)[1:end-1]

   f   = exact(0.0, 1, mesh)
   fs  = exact(0.0, 1, mesh)
   
   for n=1:nt
       
      advection_x!(mesh, f,  y, tan(0.5*dt))
      advection_y!(mesh, f, -x, sin(dt))
      advection_x!(mesh, f,  y, tan(0.5*dt))

      advection_x!(mesh, fs,  y, 0.5*dt)
      advection_y!(mesh, fs, -x, dt)
      advection_x!(mesh, fs,  y, 0.5*dt)

   end

   f, fs

end

tf, nt = 200 * π, 1000
tf, nt = 2 * pi, 10

mesh = RectMesh2D(-π, π, 256, -π, π, 256)

fe = exact(tf, nt, mesh)
@time fc, fs = with_bsl(tf, nt, mesh)
println( " errors = ", error1(fc, fe), "\t", error1(fs, fe))
x = linspace(mesh.xmin, mesh.xmax, mesh.nx+1)[1:end-1]
y = linspace(mesh.ymin, mesh.ymax, mesh.ny+1)[1:end-1]

fgnu = open("fc.dat", "w")
for i = 1:mesh.nx
   for j = 1:mesh.ny
      @printf(fgnu, "%f %f %f %f\n", x[i], y[j], fc[i,j], fe[i,j])
   end
   @printf(fgnu, "\n")
end
close(fgnu)
