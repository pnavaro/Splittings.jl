import Splittings: UniformMesh, RectMesh2D, meshgrid, poisson!

@testset "Poisson 2D on rectangular grid" begin

    mesh = RectMesh2D( UniformMesh(0, 2π, 64), UniformMesh( 0, 2π, 128))
    x = range(mesh.x1min, stop=mesh.x1max, length=mesh.n1+1)[1:end-1]
    y = range(mesh.ymin, stop=mesh.ymax, length=mesh.n2+1)[1:end-1]
    
    X, Y = meshgrid(x,y)
    
    ex   = zeros(Complex{Float64}, (mesh.n1, mesh.n2))
    ey   = zeros(Complex{Float64}, (mesh.n1, mesh.n2))
    ρ    = zeros(Complex{Float64}, (mesh.n1, mesh.n2))
    
    ρ   .= - 2 * sin.(X) .* cos.(Y);
    poisson!( ρ, mesh, ex, ey)
    @test maximum(abs.( ex .- (cos.(X) .* cos.(Y)))) ≈ 0.0 atol = 1e-15
    @test maximum(abs.( ey .+ (sin.(X) .* sin.(Y)))) ≈ 0.0 atol = 1e-15

end
