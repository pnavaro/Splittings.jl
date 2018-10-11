import Splittings: UniformMesh, RectMesh2D, meshgrid, poisson!

@testset "Poisson 2D on rectangular grid" begin

    mesh = RectMesh2D( UniformMesh(0, 2π, 64), UniformMesh( 0, 2π, 128))
    x = range(mesh.xmin, stop=mesh.xmax, length=mesh.nx+1)[1:end-1]
    y = range(mesh.ymin, stop=mesh.ymax, length=mesh.ny+1)[1:end-1]
    
    X, Y = meshgrid(x,y)
    
    ex   = zeros(Complex{Float64}, (mesh.nx, mesh.ny))
    ey   = zeros(Complex{Float64}, (mesh.nx, mesh.ny))
    ρ    = zeros(Complex{Float64}, (mesh.nx, mesh.ny))
    
    ρ   .= - 2 * sin.(X) .* cos.(Y);
    poisson!( ρ, mesh, ex, ey)
    @test maximum(abs.( ex .- (cos.(X) .* cos.(Y)))) ≈ 0.0 atol = 1e-15
    @test maximum(abs.( ey .+ (sin.(X) .* sin.(Y)))) ≈ 0.0 atol = 1e-15

end
