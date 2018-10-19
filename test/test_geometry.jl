import Splittings:Geometry

@testset " Geometry " begin

geom1 = Geometry( 11, 21, 0.0, 0.0, 0.1, 0.1 )
geom2 = Geometry( 11, 21, (0.0, 1.0), (0.0, 2.0), :none )

@test geom1.n1     == geom2.n1     == 11
@test geom1.n2     == geom2.n2     == 21
@test geom1.delta1 == geom2.delta1 == 0.1
@test geom1.delta2 == geom2.delta2 == 0.1
@test geom1.x1min  == geom2.x1min  == 0.0
@test geom1.x2min  == geom2.x2min  == 0.0
@test geom1.x1grid[end] == geom2.x1grid[end]  == 1.0
@test geom1.x2grid[end] == geom2.x2grid[end]  == 2.0

end

@testset " Geometry perx" begin

geom1 = Geometry( 10, 21, 0.0, 0.0, 0.1, 0.1 )
geom2 = Geometry( 10, 21, (0.0, 1.0), (0.0, 2.0), :perx )

@test geom1.n1     == geom2.n1     == 10
@test geom1.n2     == geom2.n2     == 21
@test geom1.delta1 == geom2.delta1 == 0.1
@test geom1.delta2 == geom2.delta2 == 0.1
@test geom1.x1min  == geom2.x1min  == 0.0
@test geom1.x2min  == geom2.x2min  == 0.0

end

@testset " Geometry pery" begin

geom1 = Geometry( 11, 20, 0.0, 0.0, 0.1, 0.1 )
geom2 = Geometry( 11, 20, (0.0, 1.0), (0.0, 2.0), :pery )

@test geom1.n1     == geom2.n1     == 11
@test geom1.n2     == geom2.n2     == 20
@test geom1.delta1 == geom2.delta1 == 0.1
@test geom1.delta2 == geom2.delta2 == 0.1
@test geom1.x1min  == geom2.x1min  == 0.0
@test geom1.x2min  == geom2.x2min  == 0.0

end

@testset " Geometry perxy" begin

geom1 = Geometry( 10, 20, 0.0, 0.0, 0.1, 0.1 )
geom2 = Geometry( 10, 20, (0.0, 1.0), (0.0, 2.0), :perxy )

@test geom1.n1     == geom2.n1     == 10
@test geom1.n2     == geom2.n2     == 20
@test geom1.delta1 == geom2.delta1 == 0.1
@test geom1.delta2 == geom2.delta2 == 0.1
@test geom1.x1min  == geom2.x1min  == 0.0
@test geom1.x2min  == geom2.x2min  == 0.0

end
