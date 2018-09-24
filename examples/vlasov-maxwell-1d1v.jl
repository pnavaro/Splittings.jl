#=

1D Vlasov–Ampere system
-----------------------
 
 $$
 \frac{\partial f}{\partial t} + \upsilon \frac{\partial f}{\partial x}
 - E(t,x) \frac{\partial f}{\partial \upsilon} = 0
 $$
 
 $$
 \frac{\partial E}{\partial t} = - J = \int f\upsilon \; d\upsilon
 $$

Algorithm 
---------
 
 - For each $j$ compute discrete Fourier transform in $x$ of 
   $f^n(x_i,\upsilon_j)$ yielding $f_k^n(\upsilon_j)$, 
 - For $ k \neq 0 $
 
     - Compute 
     
     $$f^{n+1}_k(\upsilon_j) = e^{−2i\pi k \upsilon
     \Delta t/L} f_n^k(\upsilon_j),$$
     
     - Compute 
     
     $$\rho_k^{n+1} = \Delta \upsilon \sum_j􏰄 f^{n+1}_k(\upsilon_j),$$
     
     - Compute
     
     $$E^{n+1}_k = \rho^{n+1}_k L/(2i\pi k \epsilon_0),$$
     
 - For $k = 0$ do nothing: 
 
 $$f_{n+1}(\upsilon_j) = f^n_k(\upsilon_j), E^{n+1}_k = E^n_k$$.
 
 - Perform inverse discrete Fourier transform of $E^{n+1}_k$ and for each $j$ of $f^{n+1}_k (\upsilon_j)$.

=#

using ProgressMeter, Plots, LinearAlgebra
using Splittings

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
    
    bar = Progress(nt,1)
        
    for i in 1:nt
        push!(nrj, log10(sqrt((sum(E.^2))*meshx.dx)))
        advection_v!(fᵀ, meshx, meshv, E,  0.5*dt)
        transpose!(f,fᵀ)
        E = advection_x!( f, meshx, meshv, dt)
        transpose!(fᵀ,f)
        advection_v!(fᵀ, meshx, meshv, E,  0.5*dt)
        next!(bar)
    end
    nrj
end


nx, nv = 64, 128
xmin, xmax =  0., 4*π
vmin, vmax = -6., 6.
tf = 60
nt = 600

t =  range(0,stop=tf,length=nt)
plot(t, vm1d(nx, nv, xmin, xmax, vmin, vmax, tf, nt) )
plot!(t, log10.(2.6e-3*exp.(-0.1533*t)))

@time vm1d(nx, nv, xmin, xmax, vmin, vmax, tf, nt)
