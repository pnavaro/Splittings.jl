export distibution

export GaussianProduct

struct GaussianSum

    a  :: Float64
    b  :: Float64
    xa :: Float64
    xb :: Float64
    ra :: Float64
    rb :: Float64

end

"""

    distribution( mesh1, mesh2, gauss)

    ```math
    f(x,υ) = a \exp ( \frac{x-xa}{ra}^2 ) + 
             b \exp ( \frac{υ-υb}{rb}^2 )
    ```

"""
function distribution( mesh1  :: UniformMesh,
                       mesh2  :: UniformMesh,
		       gauss  :: GaussianSum)

    x = mesh1.points
    v = mesh2.points
    a , b  = gauss.a,  gauss.b
    xa, xb = gauss.xa, gauss.xb
    ra, rb = gauss.ra, gauss.rb

    @. a * exp(-((x-xa)/ra)^2) * transpose(exp(-((v-xb)/rv)^2))

end

export GaussianProduct

struct GaussianProduct

    a  :: Float64
    ra :: Float64
    rb :: Float64

end

"""

    distribution( mesh1, mesh2, gauss)

    ```math
    f(x,υ) = a \exp (\frac{x}{ra}^2) ⋅ \exp (\frac{υ}{rb}^2)
    ```

"""
function distribution( mesh1  :: UniformMesh,
                       mesh2  :: UniformMesh,
		       gauss  :: GaussianProduct)

    x = mesh1.points
    v = mesh2.points
    rx, rv = gauss.rx, gauss.rv

    @. a * exp(-(x/rx)^2) * transpose(exp(-(v/rv)^2))

end

