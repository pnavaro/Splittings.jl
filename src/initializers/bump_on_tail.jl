struct BumpOnTail

    α  :: Float64
    kx :: Float64

end


function distribution( mesh1 :: UniformMesh,
                       mesh2 :: UniformMesh,
		       bot   :: BumpOnTail)

    x = mesh1.points
    v = mesh2.points
    α, kx = bot.α, bot.kx
    @. (1.0+α*cos(kx*x)) / (10*sqrt(2π)) * 
        transpose(9*exp(-0.5*v^2)+2*exp(-2*(v-4.5)^2))

end

