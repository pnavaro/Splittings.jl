import QuadGK:quadgk
import Roots:find_zeros


" Compute M₀ by solving F(m) = 0 "
function mag(β, mass)
    
    F(m) = begin
        g(x, n, m) = (1 / π) * (exp(β * m * cos(x)) * cos(n * x))
        bessel0(x) = g(x, 0, m) 
        bessel1(x) = g(x, 1, m)
        mass * quadgk(bessel1, 0, π)[1] / quadgk(bessel0, 0, π)[1] - m
    end
    
    find_zero(F, (0, mass))
end

#-

function Norm(f::Array{Float64,2}, dx, dv)
   return dx * sum(dv * sum(real(f), dims=1))
end

struct HMF

    b    :: Float64
    m    :: Float64
    ϵ    :: Float64

    function HMF( mass, T, ϵ )

        b = 1 / T
        m = mag(b, mass)
        w = sqrt(m)
        new( b, m, ϵ )
    end

end

function distribution( meshx :: UniformMesh,
                       meshv :: UniformMesh,
		       hmf   :: HMF)
    b  = hmf.b
    m  = hmf.m
    dx = meshx.step
    dv = meshv.step
    x  = meshx.points
    v  = meshv.points
    
    f  = transpose(exp.(-b * (v.^2 / 2))) .* exp.(b * m * cos.(x))
    a  = mass / Norm(real(f), dx, dv)
    @. a * transpose(exp(-b * (v^2) / 2)) * exp(b*m*cos(x)) * (1+ϵ*cos(x))
    
end 
