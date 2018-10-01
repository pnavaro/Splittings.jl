function compute_interpolants( n::Int, f::Array{Float64})
        
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


function interpolate( coeffs::Array{Float64,1}, nx::Int,
        xmin::Float64, xmax::Float64, x::Float64 )

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
      return (1.0/6.0)*(t2 + t4)
end
