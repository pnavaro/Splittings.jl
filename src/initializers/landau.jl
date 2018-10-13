struct Landau

    ϵ  :: Float64
    kx :: Float64

end


function distribution( meshx  :: UniformMesh,
                       meshv  :: UniformMesh,
		       landau :: Landau)

    x = meshx.points
    v = meshv.points
    ϵ, kx = landau.ϵ, landau.kx

    (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))

end

