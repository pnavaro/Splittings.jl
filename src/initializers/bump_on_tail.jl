struct BumpOnTail

    α  :: Float64
    kx :: Float64

end


function distribution( meshx :: UniformMesh,
                       meshv :: UniformMesh,
		       bot   :: BumpOnTail)

    x = meshx.points
    v = meshv.points
    α, kx = bot.α, bot.kx
    @. (1.0+α*cos(kx*x)) / (10*sqrt(2π)) * 
        transpose(9*exp(-0.5*v^2)+2*exp(-2*(v-4.5)^2))

end

