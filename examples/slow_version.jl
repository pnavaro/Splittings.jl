using FFTW
using LinearAlgebra
using Plots, ProgressMeter
using BenchmarkTools
import Splittings: Mesh, error1

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

    for n=1:nt

       for (i, xx) in enumerate(x)
           f[i,:]=real(ifft(exp.(1im*xx*ky*tan(dt/2)).*fft(f[i,:])))
       end

       for (j, yy) in enumerate(y)
           f[:,j]=real(ifft(exp.(-1im*yy*kx*sin(dt)).*fft(f[:,j])))
       end

       for (i, xx) in enumerate(x)
           f[i,:]=real(ifft(exp.(1im*xx*ky*tan(dt/2)).*fft(f[i,:])))
       end
   end

    f
end

function exact(tf, nt, mesh)

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

const tf = 200 * π
const nt = 1000

mesh = Mesh(128, 256, -π, π, -π, π)

fe = exact(tf, nt, mesh)

@time println( " error = ", error1(naive_from_matlab(tf, nt, mesh), fe))
