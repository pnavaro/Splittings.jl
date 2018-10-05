using FFTW, LinearAlgebra

"""
    bspline(p, j, x)

    Return the value at x in [0,1[ of the B-spline with
    integer nodes of degree p with support starting at j.
    Implemented recursively using the de Boor's recursion formula

"""
function bspline(p::Int, j::Int, x::Float64)
   if p == 0
       if j == 0
           return 1.0
       else
           return 0.0
       end
   else
       w = (x - j) / p
       w1 = (x - j - 1) / p
   end
   ( w       * bspline(p - 1, j    , x) +
    (1 - w1) * bspline(p - 1, j + 1, x))
end

"""
    interpolate( p, f, delta, alpha)

    Compute the interpolating spline of degree p of odd
    degree of a 1D function f on a periodic uniform mesh, at
    all points x-alpha

"""
function interpolate(p::Int, f::Vector{Float64}, delta::Float64, 
			   alpha::Float64)

   n = size(f)[1]
   modes = 2 * pi * (0:n-1) / n
   eig_bspl = zeros(Complex{Float64},n)
   eig_bspl .= bspline(p, -div(p+1,2), 0.0)
   for j in 1:div(p+1,2)-1
      eig_bspl .+= (bspline(p, j-div(p+1,2), 0.0)
         * 2 * cos.(j * modes))
   end
   ishift = floor(- alpha / delta)
   beta = - ishift - alpha / delta
   eigalpha = zeros(Complex{Float64},n)
   for j in -div(p-1,2):div(p+1,2)
      eigalpha .+= (bspline(p, j-div(p+1,2), beta)
         .* exp.((ishift + j) * 1im .* modes))
   end

   real(ifft(fft(f) .* eigalpha ./ eig_bspl))

end

"""
    advection!( mesh, f, v, dt, axis)

    Advection of a 2d function `f` discretized on a 2d `mesh`
    along the input axis at velocity `v`

"""
function advection!(f::Array{Float64,2}, mesh::RectMesh1D1V,
                    v::Any, dt::Float64; axis=0)

    @assert ( axis > 0 )

    if (axis == 1)
        for j in 1:mesh.nv
            alpha = v[j] * dt
            f[:,j] .= interpolate(3, f[:,j], mesh.dx, alpha)
        end
    else
        for i in 1:mesh.nx
            alpha = v[i] * dt
            f[i,:] .= interpolate(3, f[i,:], mesh.dv, alpha)
        end
    end

end

"""
    advection!(f, p, mesh, v, nv, dt)

    Advection of a 2d function `f` along its first dimension with
    velocity `v`. Since the fft are computed inplace, the function 
    must be represented by a Complex{Float64} 2d array.

"""
function advection!(f::Array{Complex{Float64},2}, p::Int, 
                    mesh::RectMesh1D, v::Array{Float64,1}, 
                    nv::Int, dt::Float64)

   nx = mesh.nx
   dx = mesh.dx
   modes = [2π * i / nx for i in 0:nx-1]
   # compute eigenvalues of degree p b-spline matrix
   eig_bspl = zeros(Float64, nx)
   eig_bspl .= Splittings.bspline(p, -div(p+1,2), 0.0)
   for i in 1:div(p+1,2)-1
      eig_bspl .+= Splittings.bspline(p, i - div(p+1,2), 0.0) * 2 .* cos.(i * modes)
   end
   eigalpha = zeros(Complex{Float64}, nx)

   fft!(f,1)

   for j in 1:nv
      alpha = dt * v[j] / dx

      # compute eigenvalues of cubic splines evaluated at displaced points
      ishift = floor(-alpha)
      beta   = -ishift - alpha
      fill!(eigalpha,0.0im)
      for i in -div(p-1,2):div(p+1,2)
         eigalpha .+= (Splittings.bspline(p, i-div(p+1,2), beta)
                        .* exp.((ishift+i) * 1im .* modes))
      end

      # compute interpolating spline using fft and properties of circulant matrices

      f[:,j] .*= eigalpha ./ eig_bspl

   end

   ifft!(f,1)

end

function advection!( mesh::RectMesh1D1V, f::Array{Float64,2}, dt)
    
    lx = mesh.xmax - mesh.xmin
    for j in 1:mesh.nv
        coeffs = compute_interpolants(mesh.nx, f[:,j])        
        for i in 1:mesh.nx
            x_new = mesh.x[i] - dt * mesh.v[j]
            x_new = mesh.xmin + mod(x_new - mesh.xmin, lx)
            f[i,j] = interpolate(coeffs, mesh.nx, mesh.xmin, mesh.xmax, x_new)
        end
    end
    
end

function advection!( mesh::RectMesh1D1V, f::Array{Float64,2}, e, dt)
    
    lv = mesh.vmax-mesh.vmin
    for i in 1:mesh.nx
        coeffs = compute_interpolants(mesh.nv, f[i,:])       
        for j in 1:mesh.nv
            v_new = mesh.v[j] - dt * e[i] 
            v_new = mesh.vmin + mod(v_new - mesh.vmin, lv)
            f[i,j] = interpolate(coeffs, mesh.nv, mesh.vmin, mesh.vmax, v_new)
        end
    end
    
end            
