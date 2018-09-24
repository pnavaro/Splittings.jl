using Plots, LinearAlgebra
using Splittings
pyplot()

"""

     vm1d( nx, nv, xmin, xmax, vmin, vmax , tf, nt)

Compute Landau damping by solving Vlasov-Ampere system.

# Arguments

- `nx::Integer`: the number of points along x axis.
- `nv::Integer`: the number of points along ⁠υ axis.
- `xmin::Float`: the origin of mesh along x axis.
- `xmax::Float`: the end of mesh along x axis.
- `vmin::Float`: the origin of mesh along υ axis.
- `vmax::Float`: the end of mesh along υ axis.
- `tf::Float`  : the final time of the simulation.
- `nt::Integer`: the number of time steps.

# 1D Vlasov–Ampere system

 ``
 \\frac{∂f}{∂t} + υ \\frac{∂f}{∂x}
 - E(t,x) \\frac{∂f}{∂υ} = 0
 ``

 ``
 \\frac{∂E}{∂t} = - J = ∫ fυ dυ
 ``

## Algorithm

 - For each ``j`` compute discrete Fourier transform in ``x`` of
   ``(x_i,υ_j)`` yielding ``f_k^n(υ_j)``,
 - For `` k ≂̸ 0 ``

     - Compute

     `` f^{n+1}_k(υ_j) = e^{−2iπ k υ Δt/L} f_n^k(υ_j), ``

     - Compute

     `` ρ_k^{n+1} = Δ υ ∑_j􏰄 f^{n+1}_k(υ_j), ``

     - Compute

     `` E^{n+1}_k = ρ^{n+1}_k L/(2iπkϵ_0), ``

 - For ``k = 0`` do nothing:

 `` f_{n+1}(υ_j) = f^n_k(υ_j), E^{n+1}_k = E^n_k. ``

 - Perform inverse discrete Fourier transform of ``E^{n+1}_k`` and for each
   ``j`` of ``f^{n+1}_k (υ_j)``.

# Examples

```jl
using Plots, Splittings
nx, nv = 64, 128
xmin, xmax =  0., 4π
vmin, vmax = -6., 6.
tf, nt = 60, 600
t =  range(0,stop=tf,length=nt)
plot(t, vm1d(nx, nv, xmin, xmax, vmin, vmax, tf, nt) )
plot!(t, log10.(2.6e-3*exp.(-0.1533*t)))
```

"""
function vm1d( nx, nv, xmin, xmax, vmin, vmax , tf, nt)

    meshx = UniformMesh(xmin, xmax, nx)
    meshv = UniformMesh(vmin, vmax, nv)

    # Initialize distribution function
    x = meshx.x
    v = meshv.x
    ϵ, kx = 0.001, 0.5

    f = zeros(Complex{Float64},(nx,nv))
    fᵀ= zeros(Complex{Float64},(nv,nx))

    f .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))
    transpose!(fᵀ,f)

    ρ = compute_rho(meshv, f)
    E = compute_e(meshx, ρ)

    nrj = Float64[]

    dt = tf / nt

    for i in 1:nt
        push!(nrj, log10(sqrt((sum(E.^2))*meshx.dx)))
        advection_v!(fᵀ, meshx, meshv, E,  0.5dt)
        transpose!(f,fᵀ)
        E = advection_x!( f, meshx, meshv, dt)
        transpose!(fᵀ,f)
        advection_v!(fᵀ, meshx, meshv, E,  0.5dt)
    end
    nrj
end


nx, nv = 64, 128
xmin, xmax =  0., 4π
vmin, vmax = -6., 6.
tf = 80
nt = 600

t =  range(0,stop=tf,length=nt)
plot(t, vm1d(nx, nv, xmin, xmax, vmin, vmax, tf, nt) )
#plot!(t, log10.(2.6e-3*exp.(-0.1533*t)))
plot!(t, -0.1533*t.-5.50)
