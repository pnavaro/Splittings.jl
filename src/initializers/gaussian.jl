struct GaussianSum

    a  :: Float64
    b  :: Float64
    ra :: Float64
    rb :: Float64

end

"""

    distribution( meshx, meshv, gauss)

    ```math
    f(x,υ) = a \exp ( \frac{x}{ra}^2 ) + 
             b \exp ( \frac{υ}{rb}^2 )
    ```

"""
function distribution( meshx  :: UniformMesh,
                       meshv  :: UniformMesh,
		       gauss  :: GaussianSum)

    x = meshx.points
    v = meshv.points
    a , b  = gauss.a,  gauss.b
    rx, rv = gauss.rx, gauss.rv

    @. a * exp(-(x/rx)^2) * transpose(exp(-(v/rv)^2))

end

struct GaussianProduct

    a  :: Float64
    ra :: Float64
    rb :: Float64

end

"""

    distribution( meshx, meshv, gauss)

    ```math
    f(x,υ) = a \exp (\frac{x}{ra}^2) ⋅ \exp (\frac{υ}{rb}^2)
    ```

"""
function distribution( meshx  :: UniformMesh,
                       meshv  :: UniformMesh,
		       gauss  :: GaussianProduct)

    x = meshx.points
    v = meshv.points
    rx, rv = gauss.rx, gauss.rv

    @. a * exp(-(x/rx)^2) * transpose(exp(-(v/rv)^2))

end

