struct Landau

    ϵ  :: Float64
    kx :: Float64

end


function distribution( mesh1  :: UniformMesh,
                       mesh2  :: UniformMesh,
		       landau :: Landau)

    x = mesh1.points
    v = mesh2.points
    ϵ, kx = landau.ϵ, landau.kx

    (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))

end

