#!/usr/bin/env python
# coding: utf-8

# ## 1D Vlasov–Ampere system
# 
# $$
# \frac{\partial f}{\partial t} + \upsilon \frac{\partial f}{\partial x}
# - E(t,x) \frac{\partial f}{\partial \upsilon} = 0
# $$
# 
# $$
# \frac{\partial E}{\partial t} = - J = \int f\upsilon \; d\upsilon
# $$

# ## Algorithm 
# 
# - For each $j$ compute discrete Fourier transform in $x$ of $f^n(x_i,\upsilon_j)$ yielding $f_k^n(\upsilon_j)$, 
# 
# - For $ k \neq 0 $
# 
#     - Compute 
#     
#     $$f^{n+1}_k(\upsilon_j) = e^{−2i\pi k \upsilon
#     \Delta t/L} f_n^k(\upsilon_j),$$
#     
#     - Compute 
#     
#     $$\rho_k^{n+1} = \Delta \upsilon \sum_j􏰄 f^{n+1}_k(\upsilon_j),$$
#     
#     - Compute
#     
#     $$E^{n+1}_k = \rho^{n+1}_k L/(2i\pi k \epsilon_0),$$
#     
# - For $k = 0$ do nothing: 
# 
# $$f_{n+1}(\upsilon_j) = f^n_k(\upsilon_j), E^{n+1}_k = E^n_k$$.
# 
# - Perform inverse discrete Fourier transform of $E^{n+1}_k$ and for each $j$ of $f^{n+1}_k (\upsilon_j)$.

# In[1]:


using ProgressMeter, FFTW, Plots, LinearAlgebra
using BenchmarkTools


# In[11]:


struct UniformMesh
    left   ::Float64
    right  ::Float64
    ncells ::Int32
end
delta(mesh::UniformMesh) = (mesh.right - mesh.left) / mesh.ncells
getpoints(mesh::UniformMesh) = range(mesh.left, stop=mesh.right,length=mesh.ncells+1)[1:end-1]


# In[12]:


"""
Compute charge density
ρ(x,t) = ∫ f(x,v,t) dv
"""
function compute_rho(mesh, f)    
    f̂ = fft(f,2)
    real(f̂[:,1])
end


# In[13]:


"""
Compute E using that ik*E = ρ
"""
function compute_e(mesh, ρ)
    L = mesh.right - mesh.left
    m = mesh.ncells÷2
    modes = 2π / L * [1.0;1:m-1;-m:-1]
    real(ifft(1im*fft(ρ)./modes))
end


# In[14]:


"""
Advection in υ
∂ f / ∂ t − E(x) ∂ f / ∂ υ  = 0 
"""
function advection_v!( fᵀ, meshx::UniformMesh, 
        meshv::UniformMesh, E, dt)
    
    n = meshv.ncells
    L = meshv.right - meshv.left
    k = 2π/L*[0:n÷2-1;-n÷2:-1]
    ek = exp.(-1im * dt * k * transpose(E))

    fft!(fᵀ, 1)
    fᵀ .= fᵀ .* ek
    ifft!(fᵀ, 1)
    
end


# In[15]:


function advection_x!( f, meshx::UniformMesh, 
        meshv::UniformMesh, dt)
    
    L = meshx.right - meshx.left
    m = div(meshx.ncells,2)
    k = 2π/L * [0:1:m-1;-m:1:-1]
    k̃ = 2π/L * [1;1:1:m-1;-m:1:-1]
    v = getpoints(meshv)
    ev = exp.(-1im*dt * k * transpose(v))    
    
    fft!(f,1)
    f .= f .* ev 
    Ek  = -1im * delta(meshv) * sum(f,dims=2) ./ k̃
    Ek[1] = 0.0
    ifft!(f,1)
    real(ifft(Ek))
    
end


# In[16]:


function vm1d( nx, nv, xmin, xmax, vmin, vmax , tf, nt)

    meshx = UniformMesh(xmin, xmax, nx)
    meshv = UniformMesh(vmin, vmax, nv)
            
    # Initialize distribution function
    x = getpoints(meshx)
    v = getpoints(meshv)
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
        push!(nrj, log10(sqrt((sum(E.^2))*delta(meshx))))
        advection_v!(fᵀ, meshx, meshv, E,  0.5*dt)
        transpose!(f,fᵀ)
        E = advection_x!( f, meshx, meshv, dt)
        transpose!(fᵀ,f)
        advection_v!(fᵀ, meshx, meshv, E,  0.5*dt)
        next!(bar)
    end
    nrj
end


# In[17]:


nx, nv = 64, 128
xmin, xmax =  0., 4*π
vmin, vmax = -6., 6.
tf = 60
nt = 600


# In[18]:


t =  range(0,stop=tf,length=nt)
plot(t, vm1d(nx, nv, xmin, xmax, vmin, vmax, tf, nt) )
plot!(t, log10.(2.6e-3*exp.(-0.1533*t)))


# In[19]:


@btime vm1d(nx, nv, xmin, xmax, vmin, vmax, tf, nt)


# In[ ]:




