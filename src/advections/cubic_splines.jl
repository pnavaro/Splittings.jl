"""

    compute_interpolants( n, f)

Compute interpolation coefficients

This function is a julia version of a Fortran code written by
   
Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences)
   
"""
function compute_interpolants( n::Int, f::Array{Float64})
        
    @assert (n > 27)

    num_terms = 27
    
    coeffs = zeros(Float64,n+3)
    d      = zeros(Float64,n)

    a   = sqrt((2.0+sqrt(3.0))/6.0)
    r_a = 1.0/a
    b   = sqrt((2.0-sqrt(3.0))/6.0)
    b_a = b/a

    d1 = f[1]
    coeff_tmp = 1.0
    for i in 0:num_terms-1
        coeff_tmp *= (-b_a)
        d1 += coeff_tmp*f[n-1-i]
    end

    d[1] = d1*r_a
    for i in 2:n-1
        d[i] = r_a*(f[i] - b*d[i-1])
    end
        
    d1        = d[n-1]
    coeff_tmp = 1.0
    for i in 1:num_terms
        coeff_tmp *= (-b_a)
        d1 += coeff_tmp*d[i]
    end

    coeffs[n] = d1*r_a
    
    for i in n-2:-1:1
        coeffs[i+1] = r_a*(d[i] - b*coeffs[i+2])
    end

    coeffs[1]   = coeffs[n]
    coeffs[n+1] = coeffs[2]
    coeffs[n+2] = coeffs[3]
    coeffs[n+3] = coeffs[4]
    coeffs
end

export interpolate

"""

    interpolate( coeffs, n, f)

Compute interpolatted value at `x` from coefficients get 
from `compute_interpolants`.

This function is a julia version of a Fortran code written by
   
Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences)
   
"""
function interpolate( coeffs :: Array{Float64,1}, 
                      nx     :: Int,
                      xmin   :: Float64, 
                      xmax   :: Float64, 
                      x      :: Float64 )

      rh        = (nx-1) / (xmax - xmin)
      t0        = (x-xmin)*rh
      cell      = floor(t0)
      dx        = t0 - cell
      cdx       = 1.0 - dx
      icell     = trunc(Int, cell+1)
      cim1      = coeffs[icell]
      ci        = coeffs[icell+1]
      cip1      = coeffs[icell+2]
      cip2      = coeffs[icell+3]
      t1        = 3.0*ci
      t3        = 3.0*cip1
      t2        = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
      t4        =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1

      (1.0/6.0) * (t2 + t4)

end

"""
     advection!( f, meshx, v,  dt)

Semi-lagrangian advection function of 2D distribution function represented 
by array `f`. The advection operates along `axis` (=1 is most efficient) 
with speed `v` during `dt`.

It uses cubic splines interpolation.

"""
function advection!( f::Array{Float64,2}, mesh::UniformMesh, 
                     v::Vector{Float64},  dt::Float64; axis=1)
    

    if (axis == 1) 

        @assert ( mesh.nx    == size(f)[1] )
        @assert ( size(v)[1] == size(f)[2] )
        lx = mesh.xmax - mesh.xmin
        for j in 1:size(v)[1]
            coeffs = compute_interpolants(mesh.nx, f[:,j])        
            for i in 1:mesh.nx
                x_new = mesh.x[i] - dt * v[j]
                x_new = mesh.xmin + mod(x_new - mesh.xmin, lx)
                f[i,j] = interpolate(coeffs, mesh.nx, mesh.xmin, mesh.xmax, x_new)
            end
        end

    else

        @assert ( mesh.nx    == size(f)[2] )
        @assert ( size(v)[1] == size(f)[1] )
    
        lv = meshv.xmax - meshv.xmin
        for i in 1:size(v)[1]
            coeffs = compute_interpolants(mesh.nx, f[i,:])       
            for j in 1:mesh.nx
                v_new = mesh.x[j] - dt * v[i] 
                v_new = mesh.xmin + mod(v_new - mesh.xmin, lv)
                f[i,j] = interpolate(coeffs, mesh.nx, mesh.xmin, mesh.xmax, v_new)
            end
        end

    end
    
end            
