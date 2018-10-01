using Splittings
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

import Splittings:compute_interpolants, interpolate

function interpolation_test(nx = 32)
    xmin, xmax = 0.0, 2π
    x = collect(range(xmin,stop=xmax,length=nx+1)[1:end-1])
    y = sin.(2π*x)
    x_new = zeros(Float64,nx)
    y_new = zeros(Float64,nx)
    coeffs = compute_interpolants(nx, y) # compute spline coefficients
    for (i, xi) in enumerate(x)
        x_new[i] = xi - 0.1
        x_new[i] = xmin + mod(x_new[i]-xmin,xmax-xmin) # apply periodic boundary conditions
        y_new[i] = interpolate(coeffs, nx, xmin, xmax, x_new[i])
    end
    maximum(abs.(sin.(2π*x_new) - y_new))
end

# write your own tests here
@test interpolation_test() == 0.0
