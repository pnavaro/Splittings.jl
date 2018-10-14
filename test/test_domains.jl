import Splittings:PeriodicDomain

@testset "Domains 1" begin


    x = PeriodicDomain(  0, 1, 10 )
    v = PeriodicDomain( -2, 2, 40 )

    m = x * transpose(v)

    @test m.x1min ==  0.
    @test m.x1max ==  0.9
    @test m.x2min == -2.
    @test m.x2max ==  1.9
end

@testset "Domains 2" begin

    x = PeriodicDomain( -1..1, 20 )
    v = PeriodicDomain( -2..2, 40 )

    m = transpose(x) * v

    @test m.x1min == -2.
    @test m.x1max ==  1.9
    @test m.x2min == -1.
    @test m.x2max ==  0.9

end

@testset "Domains 3" begin

    x = PeriodicDomain( -1:0.1:1.0 )
    v = PeriodicDomain( -2:0.1:2.0 )

    m = transpose(x) * v

    @test m.x1min == -2.
    @test m.x1max ==  2.0
    @test m.x2min == -1.
    @test m.x2max ==  1.0

end
