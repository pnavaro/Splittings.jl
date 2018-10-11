# ## Numerical solution of 2d rotation
# 
# $$
# \frac{d f}{dt} +  (v \frac{d f}{dx} - x \frac{d f}{dv}) = 0
# $$

" Julia function to compute exact solution "
function exact(tf, mesh::Splittings.RectMesh1D1V)

    f = zeros(Float64,(mesh.nx,mesh.nv))
    for (i, x) in enumerate(mesh.x), (j, v) in enumerate(mesh.v)
        xn=cos(tf)*x-sin(tf)*v
        vn=sin(tf)*x+cos(tf)*v
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(vn-1)*(vn-1)/0.1)
    end

    f

end

" Function to compute error "
function error1(f, f_exact)
    maximum(abs.(f .- f_exact))
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

function rotation_2d_fft(tf, nt, mesh::Splittings.RectMesh1D1V)

    dt = tf/nt

    nx = mesh.nx
    xmin, xmax = mesh.xmin, mesh.xmax
    dx = mesh.dx

    nv = mesh.nv
    vmin, vmax = mesh.vmin, mesh.vmax
    dv = mesh.dv

    kx = 2π/(xmax-xmin)*[0:nx÷2-1;nx÷2-nx:-1]
    kv = 2π/(vmax-vmin)*[0:nv÷2-1;nv÷2-nv:-1]

    f  = zeros(Complex{Float64},(nx,nv))
    f̂  = similar(f)
    fᵗ = zeros(Complex{Float64},(nv,nx))
    f̂ᵗ = similar(fᵗ)
    
    exkv = exp.( 1im * tan(dt/2) * kv .* transpose(mesh.x))
    ekxv = exp.(-1im * sin(dt)   * kx .* transpose(mesh.v))
    
    Px = plan_fft(f,  1)
    Pv = plan_fft(fᵗ, 1)
    
    f .= exact(0.0, mesh)
    
    for n=1:nt

        transpose!(fᵗ,f)
        mul!(f̂ᵗ, Pv, fᵗ)
        f̂ᵗ .= f̂ᵗ .* exkv
        ldiv!(fᵗ, Pv, f̂ᵗ)
        transpose!(f,fᵗ)
        
        mul!(f̂, Px, f)
        f̂ .= f̂ .* ekxv 
        ldiv!(f, Px, f̂)
        
        transpose!(fᵗ,f)
        mul!(f̂ᵗ, Pv, fᵗ)
        f̂ᵗ .= f̂ᵗ .* exkv
        ldiv!(fᵗ, Pv, f̂ᵗ)
        transpose!(f,fᵗ)

    end
    real(f)
end

tf, nt = 200π, 1000

mesh = Splittings.RectMesh1D1V(-π, π, 128, -π, π, 256)

fe = exact(tf, mesh)
println( " error = ", error1(rotation_2d_fft(tf, nt, mesh), fe))
