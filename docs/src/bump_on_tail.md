# Bump On Tail

```@example

using Splittings
using Plots

pyplot()


function vlasov_poisson(mesh  :: RectMesh1D1V, 
                        f     :: Array{Float64,2}, 
                        nstep :: Int64, 
                        dt    :: Float64)
    
    nrj = Float64[]
    advection!( f, mesh, 0.5dt, axis=1)
    for istep in 1:nstep
        rho = compute_rho(mesh, f)
        e   = compute_e(mesh, rho)
        advection!( f, mesh, e, dt)
        advection!( f, mesh, dt)
        push!(nrj, 0.5*log(sum(e.*e)*mesh.dx))
    end        
    nrj
    
end

α = 0.03
kx  = 0.3
xmin, xmax = 0.0, 2π / kx
nx, nv = 512, 512
vmin, vmax = -9., 9.
mesh = RectMesh1D1V(xmin, xmax, nx, vmin, vmax, nv)
f = zeros(Float64,(mesh.nx,mesh.nv))           
for (i,x) in enumerate(mesh.x), (j,v) in enumerate(mesh.v)
     f[i,j]  = (1.0+α*cos(kx*x)) / (10*sqrt(2π)) * (9*exp(-0.5*v^2)+2*exp(-2*(v-4.5)^2))
end


nstep = 500
t = range(0.0, stop=50.0, length=nstep)
dt = t[2]
@elapsed nrj = vlasov_poisson( mesh, f, nstep, dt)

plot(t, nrj, label="``\\frac{1}{2}log(∫e²dx)``")
savefig("bot-plot.png"); nothing # hide
```

![](bot-plot.png)

