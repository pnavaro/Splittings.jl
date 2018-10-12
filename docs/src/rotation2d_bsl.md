# Rotation of a gaussian distribution


```math
    \frac{df}{dt} +  (y \frac{df}{dx} - x \frac{df}{dy}) = 0
```


```@example
import Splittings: advection!, UniformMesh
import Splittings: @Strang
using Plots
pyplot(leg=false, ticks=nothing)

function with_bsl(tf::Float64, nt::Int)

   nx, ny = 64, 64
   meshx = UniformMesh(-π, π, nx)
   meshy = UniformMesh(-π, π, ny)
   x = meshx.x
   y = meshy.x

   dt = tf/nt

   f = zeros(Float64,(nx,ny))

   for (i, xp) in enumerate(x), (j, yp) in enumerate(y)
       xn = cos(tf)*xp - sin(tf)*yp
       yn = sin(tf)*xp + cos(tf)*yp
       f[i,j] = exp(-(xn-1)*(xn-1)/0.2)*exp(-(yn-1)*(yn-1)/0.2)
   end

   anim = @animate for n=1:nt
       
      @Strang(advection!( f,  meshx,  y, tan(dt), axis=1),
              advection!( f,  meshy, -x, sin(dt), axis=2))

      surface(f)
                                      
   end

   gif(anim, "rotanim.gif", fps=15)

end

f = with_bsl( 10π, 100)
```

![](rotanim.gif)
