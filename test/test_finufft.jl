"""
 Nufft object for 2d interpolation.
 It contains fft plan and 1d array to pass data to
 nufft2d subroutine from nufft package
"""
struct :: NUFFT2D

    f1d      :: Vector{Complex{Float64}}
    f1       :: Vector{Complex{Float64}}
    f2       :: Vector{Complex{Float64}}
    fcmplx   :: Array{Complex{Float64},2}
    i        :: Vector{Int64}
    j        :: Vector{Int64}
    x        :: Vector{Float64}
    y        :: Vector{Float64}
    epsnufft :: Float64
    n1       :: Int64
    n2       :: Int64
    x1min    :: Float64
    x1max    :: Float64
    x2min    :: Float64 
    x2max    :: Float64

    function NUFFT2D(n1, x1min, x1max, n2, x2min, x2max)

        x      = zeros(Float64, n1*n2)
        y      = zeros(Float64, n1*n2)
        i      = zeros(Int64,   n1*n2)
        j      = zeros(Int64,   n1*n2)
	fcmplx = zeros(Complex{Float64}, (n1,n2))
        f1     = zeros(Complex{Float64}, (n1))
        f2     = zeros(Complex{Float64}, (n2))
        f1d    = zeros(Complex{Float64}, (n1*n2))

        new( f1d, f1, f2, fcmplx, i, j, x, y, epsnufft, 
	n1, n2, x1min, x1max, x2min, x2max) 
    end 

end

"""
Compute the fft and prepare data for nufft call
Interpolate value when the fft is already computed.
"""
function interpolate!( nufft, f )

    n1 = nufft.n1
    n2 = nufft.n2
    
    fft!(fcmplx,f)
    
    m = n2÷2
    for i=1:n1
      f2=view(fcmplx,i,:)
      for j = 
      fcmplx[i,:]=f2([m+1:n2;1:m])
    end
    m = n1/2
    for j=1,n2
        f1=view(fcmplx,:,j)
        fcmplx[:,j]=f1([[(p,p=m+1,n1)], [(p,p=1,m)]])
    end
    
    xij = (x - x1min) / (x1max - x1min)
    yij = (y - x2min) / (x2max - x2min)
    
    if ( 0 < xij && xij < 1 && 0 < yij && yij < 1 )
    
      nufft.i[1] = 1
      nufft.j[1] = 1
      nufft.x[1] = xij * 6π
      nufft.y[1] = yij * 6π
    
      nufft2d2!(x, y, f1d, 1, nufft.eps, nufft.fcmplx)
    
    end

end

using FINUFFT
using Random

Random.seed!(1)

nj = 10
nk = 11
ms = 12
mt = 13
mu = 14

tol = 1e-15

# nonuniform data
x = 3π*(1.0 .- 2*rand(nj))
y = 3π*(1.0 .- 2*rand(nj))
c = rand(nj) + 1im*rand(nj)
s = rand(nk)
t = rand(nk)
u = rand(nk)
f = rand(nk) + 1im*rand(nk)

# uniform data
F2D = rand(ms, mt) + 1im*rand(ms,mt)

modevec(m) = -floor(m/2):floor((m-1)/2+1)
k1 = modevec(ms)
k2 = modevec(mt)

#int finufft2d2(BIGINT nj,FLT* xj,FLT *yj,CPX* cj,int iflag,FLT eps,
#             BIGINT ms, BIGINT mt, CPX* fk, nufft_opts opts)
#
# Type-2 2D complex nonuniform FFT.
#
#   cj[j] =  SUM   fk[k1,k2] exp(+/-i (k1 xj[j] + k2 yj[j]))      for j = 0,...,nj-1
#           k1,k2
#   where sum is over -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2,
#
# Inputs:
#   nj     number of targets (int64, aka BIGINT)
#   xj,yj     x,y locations of targets (each a size-nj FLT array) in [-3pi,3pi]
#   fk     FLT complex array of Fourier transform values (size ms*mt,
#          increasing fast in ms then slow in mt, ie Fortran ordering).
#          Along each dimension the ordering is set by opts.modeord.
#   iflag  if >=0, uses + sign in exponential, otherwise - sign (int)
#   eps    precision requested (>1e-16)
#   ms,mt  numbers of Fourier modes given in x and y (int64)
#          each may be even or odd;
#          in either case the mode range is integers lying in [-m/2, (m-1)/2].
#   opts   struct controlling options (see finufft.h)
# Outputs:
#   cj     size-nj complex FLT array of target values
#          (ie, stored as 2*nj FLTs interleaving Re, Im).
#   returned value - 0 if success, else see ../docs/usage.rst


@testset "FINUFFT 2D2" begin
    out = complex(zeros(nj))
    ref = complex(zeros(nj))
    for j=1:nj
        for ss=1:ms
            for tt=1:mt
                ref[j] += F2D[ss, tt] * exp(1im*(k1[ss]*x[j]+k2[tt]*y[j]))
            end
        end
    end
    nufft2d2!(x, y, out, 1, tol, F2D)
    relerr_2d2 = norm(vec(out)-vec(ref), Inf) / norm(vec(ref), Inf)
    @test relerr_2d2 < 1e-13
    out2 = nufft2d2(x, y, 1, tol, F2D)
    reldiff = norm(vec(out)-vec(out2), Inf) / norm(vec(out), Inf)
    @test reldiff < 1e-14            
end

