import Splittings:PeriodicDomain

@testset "Domains 1" begin


    x = PeriodicDomain(  0, 1, 10 )
    v = PeriodicDomain( -2, 2, 40 )

    m = x * transpose(v)

    @test m.xmin ==  0.
    @test m.xmax ==  0.9
    @test m.vmin == -2.
    @test m.vmax ==  1.9
end

@testset "Domains 2" begin

    x = PeriodicDomain( -1..1, 20 )
    v = PeriodicDomain( -2..2, 40 )

    m = transpose(x) * v

    @test m.xmin == -2.
    @test m.xmax ==  1.9
    @test m.vmin == -1.
    @test m.vmax ==  0.9

end

@testset "Domains 3" begin

    x = PeriodicDomain( -1:0.1:1.0 )
    v = PeriodicDomain( -2:0.1:2.0 )

    m = transpose(x) * v

    @test m.xmin == -2.
    @test m.xmax ==  2.0
    @test m.vmin == -1.
    @test m.vmax ==  1.0

end
