using LinearAlgebra, BenchmarkTools, QuadGK, Roots, FFTW
using Splittings

using Plots
pyplot()

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

function Norm(f::Array{Float64,2}, dx, dv)
   return dx * sum(dv * sum(real(f), dims=1))
end

"""
    Compute the electric hamiltonian mean field from a 
    2D distribution function
"""
function hmf_poisson!(fᵗ::Array{Float64,2}, 
        meshx::UniformMesh, 
        meshv::UniformMesh,
        ex::Array{Float64})

    nx = meshx.nx
    rho = meshv.dx .* sum(fᵗ, dims=1)[1,1:end]
    kernel = zeros(Float64, nx)
    k = π / (meshx.xmax - meshx.xmin)
    kernel[2] = k
    ex .= real(ifft(1im * fft(rho) .* kernel * 4π ))

end

"""
    Advection along first dimension using BSL 
    - To advect along dimension 2, put transpose(f) as
      input and invert meshx and meshv
"""
function bsl_advection!(meshx::UniformMesh, 
                        meshv::UniformMesh,
                        f::Array{Float64,2},
                        v::Array{Float64,1},
                        dt)
    
    @simd for j in 1:meshv.nx
        alpha = v[j] * dt
        @inbounds f[:,j] .= Splittings.interpolate(3, f[:,j], 
            meshx.dx, alpha)
    end
end

function vlasov_hmf_gauss(nbiter = 10000, deltat = 0.1)

    mass = 1.0
    T = 0.1
    meshx = UniformMesh(-π, π, 64)
    meshv = UniformMesh(-8, 8, 64)
    
    nx, dx = meshx.nx, meshx.dx
    nv, dv = meshv.nx, meshv.dx
    x, v = meshx.x, meshv.x
    X = repeat(x,1,nv)
    V = repeat(v,1,nx)'
    ϵ = 0.1
    
    b = 1 / T
    m = mag(b, mass)
   
    w   = sqrt(m)
    f   = zeros(Float64, (nx,nv))
    fᵗ  = zeros(Float64, (nv,nx))
    
    f  .= exp.(-b .* ((V.^2 / 2) - m * cos.(X)))
    a   = mass / Norm(f, dx, dv)
    @.  f =  a * exp(-b * (((V^2) / 2) - m * cos(X))) * (1 + ϵ * cos(X))
    
    ex = zeros(Float64,nx)
    hmf_poisson!(f, meshx, meshv, ex )
    T = Float64[]
    for n in 1:nbiter
    
        gamma1 = Norm(f .* cos.(X), dx, dv)
        push!(T,gamma1)
    
        """Vlasov-HMF Strang"""
        bsl_advection!(meshx, meshv, f, v, deltat / 2)
        transpose!(fᵗ, f)
        hmf_poisson!(fᵗ, meshx, meshv, ex)
        bsl_advection!(meshv, meshx, fᵗ, ex, deltat)
        transpose!(f, fᵗ)
        bsl_advection!(meshx, meshv, f, v, deltat / 2)
    
    end
    
    #Substracting from gamma its long time average
    
    Gamma1 = Norm(f .*cos.(X), dx, dv)
    T .= T .- Gamma1
    
    range(0., stop=nbiter*deltat, length=nbiter), abs.(T)
    
end 

nbiter = 1000
deltat = 0.1
@time t, T = vlasov_hmf_gauss(nbiter, deltat);

plot(t, log.(T), 
    xlabel = "t",
    ylabel = "|C[f](t)-C[f][T]|")
