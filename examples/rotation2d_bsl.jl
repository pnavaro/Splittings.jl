# # Rotation of a gaussian distribution
# 
#md # [`notebook`](https://@__NBVIEWER_ROOT_URL__notebooks/rotation_2d_bsl.ipynb)
# 
# ```math
#     \frac{df}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0
# ```
# 
# 

import Splittings: advection!, UniformMesh
import Splittings: @Strang, CubicSpline
using Plots
pyplot()

#-

function with_bsl(tf::Float64, nt::Int)

   n1, n2 = 32, 64
   mesh1 = UniformMesh(-π, π, n1)
   mesh2 = UniformMesh(-π, π, n2)
   x = mesh1.points
   y = mesh2.points

   dt = tf/nt

   f = zeros(Float64,(n1,n2))

   for (i, xp) in enumerate(x), (j, yp) in enumerate(y)
       xn = cos(tf)*xp - sin(tf)*yp
       yn = sin(tf)*xp + cos(tf)*yp
       f[i,j] = exp(-(xn-1)*(xn-1)/0.2)*exp(-(yn-1)*(yn-1)/0.2)
   end

   anim = @animate for n=1:nt
       
	   @Strang(advection!( f,  mesh1,  y, tan(dt), CubicSpline(), 1),
		   advection!( f,  mesh2, -x, sin(dt), CubicSpline(), 2)
		   )

      surface(f)
                                      
   end

   gif(anim, "rotanim.gif", fps=15); nothing #hide

end

#-

@time f = with_bsl( 2π, 20)


# ![](rotanim.gif)
