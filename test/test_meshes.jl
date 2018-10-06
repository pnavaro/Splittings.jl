import Splittings: UniformMesh

@testset "UniformMesh" begin

mesh = UniformMesh( 0, 1, 10; endpoint=false)

@test mesh.xmin == 0.0
@test mesh.dx   == 0.1
@test mesh.xmax == 1.0
@test mesh.nx   == 10

meshx = UniformMesh( 0:0.1:0.9 )

@test mesh.xmin == 0.0
@test mesh.dx   == 0.1
@test mesh.xmax == 1.0
@test mesh.nx   == 10

end
