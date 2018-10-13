# # Bump On Tail
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/bump_on_tail.ipynb),
#

import Splittings: advection!, UniformMesh
import Splittings: compute_rho, compute_e
using Plots
using LaTeXStrings

pyplot()

#-

function vlasov_poisson(meshx  :: UniformMesh, 
                        meshv  :: UniformMesh, 
                        f      :: Array{Float64,2}, 
                        nstep  :: Int64, 
                        dt     :: Float64)
    
    x = meshx.points
    v = meshv.points

    nrj = Float64[]
    for istep in 1:nstep
        advection!( f, meshx, v, 0.5dt, axis=1)
        rho = compute_rho(meshv, f)
        e   = compute_e(meshx, rho)
        advection!( f, meshv, e, dt, axis=2)
        advection!( f, meshx, v, 0.5dt, axis=1)
        push!(nrj, 0.5*log(sum(e.*e)*meshx.step))
    end        
    nrj
    
end

#-

α = 0.03
kx  = 0.3
xmin, xmax = 0.0, 2π / kx
nx, nv = 32, 64
vmin, vmax = -9., 9.
meshx = UniformMesh(xmin, xmax, nx)
meshv = UniformMesh(vmin, vmax, nv)
f = zeros(Float64,(meshx.length,meshv.length))           
for (i,x) in enumerate(meshx.points), (j,v) in enumerate(meshv.points)
    f[i,j]  = (1.0+α*cos(kx*x)) / (10*sqrt(2π)) * (9*exp(-0.5*v^2)+2*exp(-2*(v-4.5)^2))
end

#-

nstep = 500
t = range(0.0, stop=50.0, length=nstep)
dt = t[2]
@elapsed nrj = vlasov_poisson( meshx, meshv, f, nstep, dt)

#-

plot(t, nrj, label=L"\frac{1}{2} \log(∫e²dx)")
savefig("bot-plot.png"); nothing # hide

@testset "Bump On Tail" begin    #src
@test true                       #src
end                              #src

#md # ![](bot-plot.png)

