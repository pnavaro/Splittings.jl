using Plots
pyplot()

using LinearAlgebra, BenchmarkTools

masse = 1.0
T = 0.1

power(x, c, m) = c / (x^m)

"""

   Compute function to use in Bessel integral

"""
g(x, n, β, M) = (1 / π) * (exp(β * M * cos(x)) * cos(n * x))

"""

    Bessel integral

"""
Bessel(n, β, M) = sp.integrate.quad(lambda x: g(x, n, β, M), 0, π)[0]

"""

    Compute Mo by solving F(x) = 0

"""
diff(β, M) = masse * ((Bessel(1, β, M)) / (Bessel(0, β, M))) - M

"""

 Solve equation using Dichotomy

"""
function zero(diff, β, a, b) 
    if (diff(β, a) * diff(β, b) > 0)
        return 0
    end
    while (abs(a - b) > 1e-14)
        d = (a + b) / 2.
        if diff(β, d) * diff(β, a) > 0
            a = d
        else
            b = d
        end
    end
    return d
end

mag(β) = zero(diff, β, 0, masse)

function Norm(f, meshx::UniformMesh, meshv::UniformMesh)
    i1 = meshv.dx * sum(f, 1)
    i2 = meshx.dx * sum(i1)
    return i2
end


"""

    Compute the electric field from a 2D distribution function

"""
function hmf_poisson(f, meshx::UniformMesh, meshv::UniformMesh)

    nx = size(f)[1]
    @assert (nx == meshx.ncells)

    # compute rho adding neutralizing background
    rho = meshv.dx * sum(f, 2)

    N = size(rho)
    N = float(N)

    # compute Ex using that ik*Ex = rho
    modes = π * [0:nx÷2] / (meshx.xmax - meshx.xmin)
    # For the HMF case
    kernel = copy(modes)
    kernel[0:1] = 0
    kernel[2:nx // 2] = 0
    rhok = (4π / N) * rfft(rho)
    N2 = size(rhok)[1]
    N2 = float(N2)
    ex = N2 * real(irfft(1j * rhok * kernel))
    phi = N2 * real(irfft(rhok * kernel))

    return ex, phi
end


function vlasov-hmf-gauss()

    meshv = UniformMesh(-8., 8., 64)
    meshx = UniformMesh((-1) * π, π, 64)
    X = meshx.x
    V = transpose(meshv.x)
    eps = 0.1
    
    b = 1 / T
    m = mag(b)
   
    w = sqrt(m)
    eta = exp(-b * (((V^2) / 2) - m * cos(X)))
    a = masse / Norm(eta, meshx, meshv)
    f0 = a * exp(-b * (((V^2) / 2) - m * cos(X))) * (1 + eps * cos(X))
    f = f0.copy()
    ex, phi = poisson(f, meshx, meshv)
    
    nbiter = 10000
    deltat = 0.1
    T = zeros(nbiter)
    
    for n in 1:nbiter
    
        """Computation of the desired gamma for f at time n*deltat"""
        test1 = f * cos(X)
        gamma1 = Norm(test1, meshx, meshv)
        T[n] = gamma1
    
        """Vlasov-HMF Strang"""
        f = advect_x(f, meshx, meshv, deltat / 2)
        ex, phi = poisson(f, meshx, meshv)
        f = advect_v(f, ex, meshx, meshv, deltat)
        f = advect_x(f, meshx, meshv, deltat / 2)
    
    end
    
    
    """Substracting from gamma its long time average """
    Final1 = f * cos(X)
    Gamma1 = Norm(Final1, meshx, meshv)
    for n in 1:nbiter
        T[n] = T[n] - Gamma1
    end
    
    mesht = UniformMesh(deltat, (nbiter - 1) * deltat, nbiter)
    t = mesht.getpoints()
    
    """Power law fitting"""
    nd = 1001  # indice pour choix de la constante
    td = (nd - 1) * deltat  # temps correspondant
    guess = array([abs(T[nd]) * (td^3), 3.0])
    param, variance = curve_fit(power, t, T, guess, method='lm')
    Fit = param[0] / (t^param[1])
    print(param)
    print(variance)
    
    plot(t, abs(T), 'b', label="|C[f](t)-C[f][T]|")
    plot!(t, Fit, 'r', linewidth=2, label="Envelop Fitting")
    #plt.xscale('symlog')
    #plt.yscale('symlog', linthreshy=10^(-12))
    #plt.grid(True)
    #plt.xlabel("t")
    #plt.ylabel("|C[f](t)-C[f][T]|, T=nbiter*deltat")
    #plt.legend()
    #plt.title("Decay of Gamma1 for a cosine perturbation-Envelop Fitting")
    
end 
