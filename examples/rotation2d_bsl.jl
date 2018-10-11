import Splittings
using Printf

"""

    exact( tf, nt, mesh)

    Exact solution of f after rotation during time tf

"""
function exact(tf::Float64, mesh::Splittings.RectMesh1D1V)

    f = zeros(Float64,(mesh.nx,mesh.nv))
    for (i, x) in enumerate(mesh.x), (j, v) in enumerate(mesh.v)
        xn=cos(tf)*x-sin(tf)*v
        yn=sin(tf)*x+cos(tf)*v
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn-1)*(yn-1)/0.1)
    end

    f

end

"""

    error1(f, f_exact)
 
    Compute L1 error

"""
function error1(f, f_exact)
    maximum(abs.(f .- f_exact))
end


function with_bsl(tf::Float64, nt::Int, mesh::Splittings.RectMesh1D1V)

   dt = tf/nt

   f   = exact(0.0, mesh)
   fs  = exact(0.0, mesh)
   
   for n=1:nt
       
      Splittings.advection!( f,  mesh,  mesh.v, tan(0.5*dt), axis=1)
      Splittings.advection!( f,  mesh, -mesh.x, sin(dt),     axis=2)
      Splittings.advection!( f,  mesh,  mesh.v, tan(0.5*dt), axis=1)
                                      
      Splittings.advection!( fs, mesh,  mesh.v, 0.5*dt, axis=1)
      Splittings.advection!( fs, mesh, -mesh.x,     dt, axis=2)
      Splittings.advection!( fs, mesh,  mesh.v, 0.5*dt, axis=1)

   end

   f, fs

end

tf, nt = 200π, 1000

mesh = Splittings.RectMesh1D1V(-π, π, 256, -π, π, 256)

fe = exact(tf, mesh)
@time fc, fs = with_bsl(tf, nt, mesh)
println( " errors = ", error1(fc, fe), "\t", error1(fs, fe))

fgnu = open("fc.dat", "w")
for i = 1:mesh.nx
   for j = 1:mesh.ny
      @printf(fgnu, "%f %f %f %f\n", mesh.x[i], mesh.v[j], fc[i,j], fe[i,j])
   end
   @printf(fgnu, "\n")
end
close(fgnu)
