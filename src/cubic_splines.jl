"""
    This function is a julia version of a Fortran code written by
   
    Edwin Chacon Golcher
   
    Fast spline algorithm description:
   
    - data: the array whose data must be fit with the cubic spline.
    - np: (number of points, length of the data array that must be fit with
      the spline.
    - bc_type: an integer flag describing the type of boundary conditions
      desired.
    - spline_obj: the spline object to be initialized.
    This version assumes that the data are uniformly spaced.
   
    The idea behind this spline implementation is the factorization of the
    array:
   
             + 4  1              1 +
             | 1  4  1             |
             |    .  .  .          |
             |       .  .  .       |
             |          .  .  .    |
             |             1  4  1 |
             + 1              1  4 +
   
    in the form:
   
               A = L*L^t
   
    where:
   
             + a                 b +
             | b  a                |
             |    .  .             |
      L =    |       .  .          |    (zeros not shown)
             |          .  .       |
             |             b  a    |
             +                b  a +
   
    This factorization is achieved for the values:
    a² = (2+sqrt(3))/6, b² = (2-sqrt(3))/6.
   
    The spline coefficients C are thus the solution of:
   
                     L*L^t*C = F.
   
    Hence, we solve first for D in:
   
                     L*D = F
   
    This solution is achieved in a two-step process. First we compute the
    first term d1, which is given by:
   
    d1 =
   
    1           b           b  2                   i  b  i
   ---(f(x1) - ---f(xN) + (---)*f(xN-1) + ... + (-1)(---) (f(xN-i+1)-bd(N-i)))
    a           a           a                         a
   
    The series converges (since a > b) and we can approximate the series
    by a partial sum:
   
          1                            b
    d1 = ---(f(x1) + SUM(i=1,M)(-1)^i(---)^i*f(xN-i+1))
          a                            a
   
    The rest of the terms can be found by:
   
    d(i) = 1/a*(f(x i) - b*d(i-1))
   
    Once D is known, the same procedure can be used to compute C in
   
                      L^t*C = D
   
    The same process is carried out. First, the last coefficient is
    calculated by:
   
         c(N) = 1/a*(d(N) + SUM(i=1,M) (-1)^i*(b/a)^i*d(i))
   
    And the rest of the terms, starting with C(N-1) and working backwards:
   
     c(i) = 1/a*(d(i) - b*c(i+1))

    The algorithm above is not implemented whenever the number of points is
    smaller than 28. In such cases we fall back to a more straightforward but
    also more costly implementation using a standard tridiagonal system
    solver.
   
    In the periodic case, we are solving the problem A*C=F, where
   
                   + 4  1              1 +         1
                   | 1  4  1             |
                1  |    .  .  .          |         .
          A =  --- |       .  .  .       |         .
                6  |          .  .  .    |         .
                   |             1  4  1 |
                   + 1              1  4 +      num_points
   
    In the Hermite BC case, the problem is modified as follows:
   
       + 4  2                + + c_0    +   + 6*f_0+2*delta*f_0´          +  0
       | 1  4  1             | |        |   | 6*f_1                       |  1
       |    .  .  .          | |    .   |   |    .                        |  .
    A =|       .  .  .       |*|    .   | = |    .                        |  .
       |          .  .  .    | |    .   |   |    .                        |  .
       |             1  4  1 | |        |   | 6*f_np                      |  np
       +                2  4 + +c_(np+1)+   + 6*f_(np+1)+2*delta*f_(np+1)'+ np+1

"""
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
