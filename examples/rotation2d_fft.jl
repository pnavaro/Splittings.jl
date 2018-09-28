# ## Numerical solution of 2d rotation
# 
# $$
# \frac{d f}{dt} +  (v \frac{d f}{dx} - x \frac{d f}{dv}) = 0
# $$

using  FFTW, LinearAlgebra, BenchmarkTools

import Splittings:Mesh

" Julia function to compute exact solution "
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

" Function to compute error "
function error1(f, f_exact)
    maximum(abs.(f - f_exact))
end

import FFTW: FFTWPlan

struct Advector

    f  :: Array{Complex{Float64},2}
    p1 :: FFTWPlan
    p2 :: FFTWPlan

    function Advector( f  :: Array{Complex{Float64},2},
	               fᵗ :: Array{Complex{Float64},2})

        p1 = plan_fft(f,  1)
        p2 = plan_fft(fᵗ, 1)

    end
end

function rotation_2d_fft(tf, nt, mesh::Mesh)

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

const tf = 200 * π
const nt = 1000

mesh = Mesh(128, 256, -π, π, -π, π)

fe = exact(tf, nt, mesh);

@time rotation_2d_fft(tf, nt, mesh)
@time rotation_2d_fft(tf, nt, mesh)
@time rotation_2d_fft(tf, nt, mesh)
@time rotation_2d_fft(tf, nt, mesh)
println( " error = ", error1(rotation_2d_fft(tf, nt, mesh), fe))
