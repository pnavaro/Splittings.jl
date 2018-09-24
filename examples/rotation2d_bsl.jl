using Plots, Splittings

"""

    exact( tf, nt, mesh)

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

    error1(f, f_exact)
 
    Compute L1 error

"""
function error1(f, f_exact)
    maximum(abs.(f - f_exact))
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

tf, nt = 200π, 1000

mesh = RectMesh2D(-π, π, 256, -π, π, 256)

fe = exact(tf, nt, mesh)
@time fc, fs = with_bsl(tf, nt, mesh)
println( " errors = ", error1(fc, fe), "\t", error1(fs, fe))
x = range(mesh.xmin, stop=mesh.xmax, length=mesh.nx+1)[1:end-1]
y = range(mesh.ymin, stop=mesh.ymax, length=mesh.ny+1)[1:end-1]

fgnu = open("fc.dat", "w")
for i = 1:mesh.nx
   for j = 1:mesh.ny
      @printf(fgnu, "%f %f %f %f\n", x[i], y[j], fc[i,j], fe[i,j])
   end
   @printf(fgnu, "\n")
end
close(fgnu)
