import Splittings: UniformMesh

@testset "UniformMesh" begin

mesh = UniformMesh( 0, 1, 11; endpoint=true)

@test mesh.start   == 0.0
@test mesh.step    == 0.1
@test mesh.stop    == 1.0
@test mesh.length  == 11

mesh = UniformMesh( 0, 1, 10; endpoint=false)

@test mesh.start   == 0.0
@test mesh.step    == 0.1
@test mesh.stop    == 1.0
@test mesh.length  == 10

end
