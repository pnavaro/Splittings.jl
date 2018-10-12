# Rotation of a gaussian distribution


```math
    \frac{df}{dt} +  (y \frac{df}{dx} - x \frac{df}{dy}) = 0
```


```@example
using Splittings

" Exact solution of f after rotation during time tf "
function exact(tf::Float64, mesh::RectMesh1D1V)

    f = zeros(Float64,(mesh.nx,mesh.nv))
    for (i, x) in enumerate(mesh.x), (j, v) in enumerate(mesh.v)
        xn=cos(tf)*x-sin(tf)*v
        yn=sin(tf)*x+cos(tf)*v
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn-1)*(yn-1)/0.1)
    end

    f

end

function with_bsl(tf::Float64, nt::Int, mesh::RectMesh1D1V)

   dt = tf/nt

   f  = exact(0.0, mesh)
   
   for n=1:nt
       
      @Strang(advection!( f,  mesh,  mesh.v, tan(dt), axis=1),
              advection!( f,  mesh, -mesh.x, sin(dt), axis=2))
                                      
   end

   f

end

tf, nt = 10π, 100

mesh = RectMesh1D1V(-π, π, 64, -π, π, 64)

@time f = with_bsl(tf, nt, mesh)
```
