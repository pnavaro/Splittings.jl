import Splittings:UniformMesh


"""
   Julia function to compute exact solution "

`` \\frac{d f}{dt} +  (y \\frac{df}{dx} - x \\frac{df}{dy}) = 0 ``

"""
function exact(tf, meshx::UniformMesh, meshy::UniformMesh)

    f = zeros(Float64,(meshx.nx,meshy.nx))
    for (i, x) in enumerate(meshx.x), (j, y) in enumerate(meshy.x)
        xn = cos(tf) * x - sin(tf) * y
        yn = sin(tf) * x + cos(tf) * y
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn-1)*(yn-1)/0.1)
    end

    f

end

" Function to compute error "
function error1(f, f_exact)
    maximum(abs.(f .- f_exact))
end

import FFTW: FFTWPlan, plan_fft, fft!, ifft!

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

function rotation_2d_fft(tf, nt, meshx::UniformMesh, meshy::UniformMesh)

    dt = tf/nt

    nx = meshx.nx
    xmin, xmax = meshx.xmin, meshx.xmax
    dx = meshx.dx

    ny = meshy.nx
    ymin, ymax = meshy.xmin, meshy.xmax
    dy = meshy.dx

    kx = 2π/(xmax-xmin)*[0:nx÷2-1;nx÷2-nx:-1]
    ky = 2π/(ymax-ymin)*[0:ny÷2-1;ny÷2-ny:-1]

    f  = zeros(Complex{Float64},(nx,ny))
    f̂  = similar(f)
    fᵗ = zeros(Complex{Float64},(ny,nx))
    f̂ᵗ = similar(fᵗ)
    
    exky = exp.( 1im * tan(dt/2) * ky .* transpose(meshx.x))
    ekxy = exp.(-1im * sin(dt)   * kx .* transpose(meshy.x))
    
    Px = plan_fft(f,  1)
    Py = plan_fft(fᵗ, 1)
    
    f .= exact(0.0, meshx, meshy)
    
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

tf, nt = 200π, 1000

meshx = UniformMesh(-π, π, 128; endpoint=false)
meshy = UniformMesh(-π, π, 256; endpoint=false)

fc = rotation_2d_fft(tf, nt, meshx, meshy)
fe = exact(tf, meshx, meshy)

@testset "Rotation test with Fourier advections " begin

println(error1(fc, fe))
@test rotation_2d_fft(tf, nt, meshx, meshy) ≈ fe atol = 1e-10

end
