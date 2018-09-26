# ## Numerical solution of 2d rotation
# 
# $$
# \frac{d f}{dt} +  (v \frac{d f}{dx} - x \frac{d f}{dv}) = 0
# $$

using  FFTW
using  LinearAlgebra
using  Plots, ProgressMeter
using  BenchmarkTools

# ## Original matlab code
# 
# This code is not optimized and written quickly but it runs during less than 3 seconds

# ```matlab
# tf=200*pi;Nt=1000;dt=tf/Nt;
# Nx=64;Ny=128;
# xmin=-pi; xmax=pi; dx=(xmax-xmin)/Nx;
# x=xmin:dx:xmax-dx;
# ymin=-pi; ymax=pi; dy=(ymax-ymin)/Ny;
# y=ymin:dy:ymax-dy;
# 
# kx=2*pi/(xmax-xmin)*[0:Nx/2-1,Nx/2-Nx:Nx-1-Nx];
# ky=2*pi/(ymax-ymin)*[0:Ny/2-1,Ny/2-Ny:Ny-1-Ny];
# 
# % initial condition 
# f=zeros(Nx,Ny);
# for i=1:Nx
#     xx=xmin+(i-1)*dx;
#     for j=1:Ny
#         yy=ymin+(j-1)*dy;
#         f(i,j)=exp(-(xx-1)*(xx-1)/0.1)*exp(-(yy-1)*(yy-1)/0.1);
#     end
# end
# 
# fnx=zeros(1,Nx);ffx=zeros(1,Nx);fny=zeros(1,Ny);ffy=zeros(1,Ny);
# 
# for n=1:Nt
#
#     for i=1:Nx
#         xx=xmin+(i-1)*dx;
#         ffy=fft(f(i,:));
#         fny=real(ifft(exp(sqrt(-1)*xx*ky*tan(dt/2)).*ffy));
#         f(i,:)=fny;
#     end
#     
#     for j=1:Ny
#         yy=ymin+(j-1)*dy;
#         ffx=fft(f(:,j));
#         fnx=real(ifft(exp(-sqrt(-1)*yy*kx*sin(dt)).*transpose(ffx)));
#         f(:,j)=fnx;
#     end
# 
#     for i=1:Nx
#         xx=xmin+(i-1)*dx;
#         ffy=fft(f(i,:));
#         fny=real(ifft(exp(sqrt(-1)*xx*ky*tan(dt/2)).*ffy));
#         f(i,:)=fny;
#     end
#     
# end
# 
# % compute exact solution
# fexact=zeros(Nx,Ny);
# for i=1:Nx
#     xx=xmin+(i-1)*dx;    
#     for j=1:Ny
#         yy=ymin+(j-1)*dy;
#         xn=cos(tf)*xx-sin(tf)*yy;
#         yn=sin(tf)*xx+cos(tf)*yy;
#         f_exact(i,j)=exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn-1)*(yn-1)/0.1);
#     end
# end
# % compute errors in Linfty norm
# error1=max(max(f-f_exact))
# 
# ```

# ### Julia type for mesh information

import Splittings:Mesh

# ### Julia function to compute exact solution

function exact(tf, nt, mesh::Mesh)

    dt = tf/nt
    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    dx = (xmax - xmin) / nx
    x = range(xmin, stop=xmax-dx, length=nx )

    ny = mesh.ny
    ymin, ymax = mesh.ymin, mesh.ymax
    dy = (ymax -ymin) / ny
    y  = range(ymin, stop=ymax-dy, length=ny)

    f = zeros(Float64,(nx,ny))
    for (i, xx) in enumerate(x), (j, yy) in enumerate(y)
        xn=cos(tf)*xx-sin(tf)*yy
        yn=sin(tf)*xx+cos(tf)*yy
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn-1)*(yn-1)/0.1)
    end

    f
end


#=
## Create the gif to show what we are computing

function create_gif()
    x = y = range(-π, stop=π, length=40)
    n = 100
    
    # create a progress bar for tracking the animation generation
    prog = Progress(n,1)
    
    @gif for t in range(0, stop=2π, length=n)
        f(x,y) = exp(-((cos(t)*x-sin(t)*y)-1)^2/0.2)*exp(-((sin(t)*x+cos(t)*y)-1)^2/0.2)
        
        p = plot(x, y, f, st = [:surface])
    
        plot!(p[1])
        plot!(zlims=(-0.01,1.01))
    
        next!(prog) # increment the progress bar
    end
end

create_gif()
=#

# ## Function to compute error 
# 
# - In julia the max value of an array is `maximum`.

function error1(f, f_exact)
    maximum(abs.(f - f_exact))
end


# ### Naive translation from the matlab code

function naive_from_matlab(tf, nt, mesh::Mesh)

    dt = tf/nt
    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    dx = (xmax - xmin) / nx
    x = range(xmin, stop=xmax-dx, length=nx)

    ny = mesh.ny
    ymin, ymax = mesh.ymin, mesh.ymax
    dy = (ymax -ymin) / ny
    y  = range(ymin, stop=ymax-dy, length=ny)

    kx = 2π/(xmax-xmin)*[0:nx/2-1;nx/2-nx:-1]
    ky = 2π/(ymax-ymin)*[0:ny/2-1;ny/2-ny:-1]

    f = exact(0.0, 1, mesh)

    fnx = zeros(Complex{Float64},(1,nx))
    ffx = zeros(Complex{Float64},(1,nx))
    fny = zeros(Complex{Float64},(1,ny))
    ffy = zeros(Complex{Float64},(1,ny))

    for n=1:nt
       
       for (i, xx) in enumerate(x)
           ffy    .= fft(f[i,:])
           fny    .= real(ifft(exp.(1im*xx*ky*tan(dt/2)) .* ffy))
           f[i,:] .= fny
       end
       
       for (j, yy) in enumerate(y)
           ffx    .= fft(f[:,j])
           fnx    .= real(ifft(exp.(-1im*yy*kx*sin(dt)) .* ffx'))
           f[:,j] .= fnx
       end
       
       for (i, xx) in enumerate(x)
           ffy    .= fft(f[i,:])
           fny    .= real(ifft(exp.(1im*xx*ky*tan(dt/2)) .* ffy))
           f[i,:] .= fny
       end
   end

    f
end


# ###  Vectorized version
# 
# - We remove the for loops over direction x and y by creating the 2d arrays `exky` and `ekxy`.
# - We save cpu time by computing them before the loop over time

function vectorized(tf, nt, mesh::Mesh)

    dt = tf/nt
    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    dx = (xmax - xmin) / nx
    x = range(xmin, stop=xmax-dx, length=nx)

    ny = mesh.ny
    ymin, ymax = mesh.ymin, mesh.ymax
    dy = (ymax -ymin) / ny
    y  = range(ymin, stop=ymax-dy, length=ny)

    kx = 2*π/(xmax-xmin)*[0:nx÷2-1;nx÷2-nx:-1]
    ky = 2*π/(ymax-ymin)*[0:ny÷2-1;ny÷2-ny:-1]

    f = zeros(Float64,(nx,ny))
    f .= exact(0.0, 1, mesh)

    exky = exp.(1im*tan(dt/2) * x * transpose(ky))
    ekxy = exp.(-1im*sin(dt) * kx * transpose(y))
    
    for n=1:nt
        f = real(ifft(exky .* fft(f, 2), 2))
        f = real(ifft(ekxy .* fft(f, 1), 1))
        f = real(ifft(exky .* fft(f, 2), 2))
    end

    f
end


# ## Inplace computation 
# - We remove the Float64-Complex128 conversion by allocating the distribution function `f` as a Complex array
# - Note that we need to use the inplace assignement operator ".="  to initialize the `f` array.
# - We use inplace computation for fft with the "bang" operator `!`

function inplace(tf, nt, mesh::Mesh)

    dt = tf/nt
    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    dx = (xmax - xmin) / nx
    x = range(xmin, stop=xmax-dx, length=nx)

    ny = mesh.ny
    ymin, ymax = mesh.ymin, mesh.ymax
    dy = (ymax -ymin) / ny
    y  = range(ymin, stop=ymax-dy, length=ny)

    kx = 2π/(xmax-xmin)*[0:nx÷2-1;nx÷2-nx:-1]
    ky = 2π/(ymax-ymin)*[0:ny÷2-1;ny÷2-ny:-1]

    f = zeros(Complex{Float64},(nx,ny))
    f .= exact(0.0, 1, mesh)

    exky = exp.(1im*tan(dt/2) * x * transpose(ky))
    ekxy = exp.(-1im*sin(dt) * kx * transpose(y))
    
    for n=1:nt
        fft!(f, 2)
        f .= exky .* f
        ifft!(f,2)
        fft!(f, 1)
        f .= ekxy .* f
        ifft!(f, 1)
        fft!(f, 2)
        f .= exky .* f
        ifft!(f,2)        
    end

    real(f)
end

# ### Use plans for fft

# - When you apply multiple fft on array with same shape and size, it is recommended to use fftw plan to improve computations.
# - Let's try to initialize our two fft along x and y with plans.
# - We set back `f` to Array{Float64,2}
# - We use plan_rfft and plan_irfft, they need some preallocations.

function with_fft_plans(tf, nt, mesh::Mesh)

    dt = tf/nt
    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    dx = (xmax - xmin) / nx
    x = range(xmin, stop=xmax-dx, length=nx)

    ny = mesh.ny
    ymin, ymax = mesh.ymin, mesh.ymax
    dy = (ymax -ymin) / ny
    y  = range(ymin, stop=ymax-dy, length=ny)

    kx = 2π/(xmax-xmin)*collect(0.0:1.0:nx÷2)
    ky = 2π/(ymax-ymin)*collect(0.0:1.0:ny÷2)
   
    f = zeros(Float64,(nx,ny))
    f̂̂₁= zeros(Complex{Float64},(div(nx,2)+1,ny))
    f̂₂= zeros(Complex{Float64},(nx,div(ny,2)+1))
    
    exky = exp.(1im*tan(dt/2) * x * transpose(ky))
    ekxy = exp.(-1im*sin(dt) * kx * transpose(y))
    
    fw_x = plan_rfft(f, 1)
    bw_x = plan_irfft(f̂̂₁,nx,  1)
    
    fw_y = plan_rfft(f, 2)
    bw_y = plan_irfft(f̂₂,ny,  2)
    
    f .= exact(0.0, 1, mesh)
    
    for n=1:nt
        
        f̂₂  = fw_y * f
        f̂₂ .= f̂₂ .* exky
        f   = bw_y * f̂₂
        
        f̂₁  = fw_x * f
        f̂₁ .= f̂₁ .* ekxy 
        f   = bw_x * f̂₁
        
        f̂₂  = fw_y * f
        f̂₂ .= f̂₂.* exky
        f   = bw_y * f̂₂
        
    end

    f
end

# ## Inplace computation and fft plans
# 
# - In this example we use both plans and inplace computation
# 
# To apply fft plan to an array A, we use a preallocated output array Â by calling `mul!(Â, plan, A)`. 
# The input array A must be a complex floating-point array like the output Â.
# The inverse-transform is computed inplace by applying `inv(P)` with `ldiv!(A, P, Â)`.

function with_fft_plans_inplace(tf, nt, mesh::Mesh)

    dt = tf/nt
    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    dx = (xmax - xmin) / nx
    x = range(xmin, stop=xmax-dx, length=nx)

    ny = mesh.ny
    ymin, ymax = mesh.ymin, mesh.ymax
    dy = (ymax -ymin) / ny
    y  = range(ymin, stop=ymax-dy, length=ny)

    kx = 2π/(xmax-xmin)*[0:nx/2-1;nx/2-nx:-1]
    ky = 2π/(ymax-ymin)*[0:ny/2-1;ny/2-ny:-1]

    f = zeros(Complex{Float64},(nx,ny))
    f̂ = zeros(Complex{Float64},(nx,ny))
    
    exky = exp.( 1im * tan(dt/2) * x  * transpose(ky))
    ekxy = exp.(-1im * sin(dt)   * kx * transpose(y))
    
    Px = plan_fft(f, 1)    
    Py = plan_fft(f, 2)
    
    f .= exact(0.0, 1, mesh)
    
    for n=1:nt
        
        mul!(f̂, Py, f)
        f̂ .= f̂ .* exky
        ldiv!(f, Py, f̂)
        
        mul!(f̂, Px, f)
        f̂ .= f̂ .* ekxy 
        ldiv!(f, Px, f̂)
        
        mul!(f̂, Py, f)
        f̂ .= f̂ .* exky
        ldiv!(f, Py, f̂)
        
    end

    real(f)
end

# ## Explicit transpose of `f`
# 
# - Multidimensional arrays in Julia are stored in column-major order.
# - FFTs along y are slower than FFTs along x
# - We can speed-up the computation by allocating the transposed `f` 
# and transpose f for each advection along y.

function with_fft_transposed(tf, nt, mesh::Mesh)

    dt = tf/nt
    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    dx = (xmax - xmin) / nx
    x = range(xmin, stop=xmax-dx, length=nx)

    ny = mesh.ny
    ymin, ymax = mesh.ymin, mesh.ymax
    dy = (ymax -ymin) / ny
    y  = range(ymin, stop=ymax-dy, length=ny)

    kx = 2π/(xmax-xmin)*[0:nx/2-1;nx/2-nx:-1]
    ky = 2π/(ymax-ymin)*[0:ny/2-1;ny/2-ny:-1]

    f  = zeros(Complex{Float64},(nx,ny))
    f̂  = similar(f)
    fᵗ = zeros(Complex{Float64},(ny,nx))
    f̂ᵗ = similar(fᵗ)
    
    exky = exp.( 1im * tan(dt/2) * ky * transpose(x))
    ekxy = exp.(-1im * sin(dt)   * kx * transpose(y))
    
    Px = plan_fft(f,  1)    
    Py = plan_fft(fᵗ, 1)
    
    f .= exact(0.0, 1, mesh)
    
    for n=1:nt
        transpose!(fᵗ,f)
        mul!(f̂ᵗ, Py, fᵗ)
        f̂ᵗ .= f̂ᵗ .* exky
        ldiv!(fᵗ, Py, f̂ᵗ)
        transpose!(f,fᵗ)
        
        mul!(f̂, Px, f)
        f̂ .= f̂ .* ekxy 
        ldiv!(f, Px, f̂)
        
        transpose!(fᵗ,f)
        mul!(f̂ᵗ, Py, fᵗ)
        f̂ᵗ .= f̂ᵗ .* exky
        ldiv!(fᵗ, Py, f̂ᵗ)
        transpose!(f,fᵗ)
    end
    real(f)
end

d = Dict() 
const tf = 200 * π
const nt = 1000

mesh = Mesh(128, 256, -π, π, -π, π)

fe = exact(tf, nt, mesh);

@time println( " error = ", error1(naive_from_matlab(tf, nt, mesh), fe))
@time println( " inplace = ", error1(inplace(tf, nt, mesh), fe))
@time println( " vectorized = ", error1(vectorized(tf, nt, mesh), fe))
@time println( " with_fft_plans = ", error1(with_fft_plans(tf, nt, mesh), fe))
@time println( " with_fft_plans_inplace = ", error1(with_fft_plans_inplace(tf, nt, mesh), fe))
@time println( " with_fft_transposed = ", error1(with_fft_transposed(tf, nt, mesh), fe))

inplace_bench = @benchmark inplace(tf, nt, mesh)
vectorized_bench = @benchmark vectorized(tf, nt, mesh)
with_fft_plans_bench = @benchmark with_fft_plans(tf, nt, mesh)
with_fft_plans_inplace_bench = @benchmark with_fft_plans_inplace(tf, nt, mesh)
with_fft_transposed_bench = @benchmark with_fft_transposed(tf, nt, mesh)

d["vectorized"] = minimum(vectorized_bench.times) / 1e6
d["inplace"] = minimum(inplace_bench.times) / 1e6
d["with_fft_plans"] = minimum(with_fft_plans_bench.times) / 1e6
d["with_fft_plans_inplace"] = minimum(with_fft_plans_inplace_bench.times) / 1e6
d["with_fft_transposed"] = minimum(with_fft_transposed_bench.times) / 1e6;

for (key, value) in sort(collect(d), by=last)
    println(rpad(key, 25, "."), lpad(round(value, 1), 6, "."))
end

# ## Conclusion
# - Use pre-allocations of memory and inplace computation are very important
# - Try to always do computation on data contiguous in memory
# - In this notebook, time consumed in compilation time is not measured. 
# If we consider bigger and longer problem, it is not an issue.
