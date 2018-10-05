using Test
using IntervalSets

import Splittings:compute_interpolants, interpolate

@testset "CubicSpline interpolation" begin

    function interpolation_test(nx = 128)
    
        xmin, xmax = 0.0, 1.0
        x = collect(range(xmin,stop=xmax,length=nx))
        y = sin.(2π*x)
        x_new = zeros(Float64,nx)
        y_new = zeros(Float64,nx)
        coeffs = compute_interpolants(nx, y) 
        for (i, xi) in enumerate(x)
            x_new[i] = xi - 0.1
            x_new[i] = xmin + mod(x_new[i]-xmin,xmax-xmin) 
            y_new[i] = interpolate(coeffs, nx, xmin, xmax, x_new[i])
        end
        maximum(abs.(sin.(2π*(x.-0.1)) - y_new))
    
    end

    @test ≈(interpolation_test(), 0.0, atol=1e-7)

end

@testset "Domains 1" begin

    import Splittings:PeriodicDomain

    x = PeriodicDomain( -1, 1, 21 )
    v = PeriodicDomain( -2, 2, 41 )

    m = x * transpose(v)

    @test m.xmin == -1.
    @test m.xmax ==  0.9
    @test m.vmin == -2.
    @test m.vmax ==  1.9
end

@testset "Domains 2" begin

    x = PeriodicDomain( -1..1, 21 )
    v = PeriodicDomain( -2..2, 41 )

    m = transpose(x) * v

    @test m.xmin == -2.
    @test m.xmax ==  1.9
    @test m.vmin == -1.
    @test m.vmax ==  0.9

end
