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

@testset "2D2" begin
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

