import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
from scipy import integrate
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.pyplot import pause
from scipy.optimize import curve_fit
import time

import sys

masse = 1.0
T = 0.1

def power(x, c, m):
    return c / (x ** m)

def g(x, n, beta, M):  # calcul de la fonction a integrer dans l'integrale Bessel
    return (1 / np.pi) * (np.exp(beta * M * np.cos(x)) * np.cos(n * x))

def Bessel(n, beta, M):  # integrale de Bessel
    return sp.integrate.quad(lambda x: g(x, n, beta, M), 0, np.pi)[0]

def diff(beta, M):  # Pour resoudre l'equation sur Mo sous la forme F(x) = 0
    return masse * ((Bessel(1, beta, M)) / (Bessel(0, beta, M))) - M

def zero(diff, beta, a, b):  # On resoud l'equation par dichotomie
    if diff(beta, a) * diff(beta, b) > 0:
        # print('pas de M_{0} solution entre '+str(a)+' et '+str(b)+'!')
        return 0
    while abs(a - b) > 1e-14:
        d = (a + b) / 2.
        if diff(beta, d) * diff(beta, a) > 0:
            a = d
        else:
            b = d
    # print('la solution M_{0} est '+str(d))
    return d

def mag(beta):
    return zero(diff, beta, 0, masse)

def Norm(f, meshx, meshv):
    i1 = meshv.deltax * np.sum(f, axis=0)
    i2 = meshx.deltax * np.sum(i1)
    return i2

class UniformMesh:
    def __init__(self, xmin, xmax, ncells):
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.deltax = (self.xmax - self.xmin) / ncells
        self.ncells = ncells

    def getx(self, i):
        return self.xmin + i * self.deltax

    def getpoints(self):
        return self.xmin + self.deltax * np.arange(self.ncells)

def bspline(p, j, x):
    """Return the value at x in [0,1[ of the B-spline with integer nodes of degree p with support starting at j.
    Implemented recursively using the de Boor's recursion formula"""
    x = float(x)
    assert ((x >= 0.0) & (x <= 1.0))
    assert ((type(p) == int) & (type(j) == int))
    if p == 0:
        if j == 0:
            return 1.0
        else:
            return 0.0
    else:
        w = (x - j) / p
        w1 = (x - j - 1) / p
        return w * bspline(p - 1, j, x) + (1 - w1) * bspline(p - 1, j + 1, x)

def interpSpline(p, f, mesh, alpha):
    """compute the interpolating spline of degree p of odd degree of a function f on a periodic uniform mesh, at
    all points xi-alpha"""
    assert p & 1  # check that p is odd
    n = np.size(f)
    assert (n == mesh.ncells)
    # compute eigenvalues of degree p b-spline matrix
    modes = 2 * np.pi * (np.arange(n)) / n
    eig_bspl = bspline(p, -(p + 1) // 2, 0.0)
    for j in range(1, (p + 1) // 2):
        eig_bspl += bspline(p, j - (p + 1) // 2, 0.0) * 2 * np.cos(j * modes)
    # compute eigenvalues of cubic splines evaluated at displaced points
    ishift = np.floor(-alpha / mesh.deltax)
    beta = -ishift - alpha / mesh.deltax
    eigalpha = np.zeros(n, dtype=complex)
    for j in range(-(p - 1) // 2, (p + 1) // 2 + 1):
        eigalpha += bspline(p, j - (p + 1) // 2, beta) * np.exp((ishift + j) * 1j * modes)
    # compute interpolating spline using fft and properties of circulant matrices
    interpSpline = np.real(np.fft.ifft(np.fft.fft(f) * eigalpha / eig_bspl))
    return interpSpline


def poisson(f, meshx, meshv):
    """compute the electric field from a 2D distribution function"""
    nx = np.size(f, 1)
    assert (nx == meshx.ncells)
    # compute rho adding neutralizing background
    rho = meshv.deltax * np.sum(f, axis=0)

    N = np.size(rho)
    N = float(N)

    # compute Ex using that ik*Ex = rho
    modes = np.pi * (np.arange(nx // 2 + 1)) / (meshx.xmax - meshx.xmin)
    # For the HMF case
    kernel = modes.copy()
    kernel[0:1] = 0
    kernel[2:nx // 2] = 0
    rhok = (4 * np.pi / N) * np.fft.rfft(rho)
    N2 = np.size(rhok)
    N2 = float(N2)
    # For Vlasov Poisson
    # modes[0] = 1.  # avoid division by 0
    # phi = real(fft.irfft(rhok/(modes**2)))
    # ex = real(fft.irfft(-1j*rhok/modes))
    # For HMF
    ex = N2 * np.real(np.fft.irfft(1j * rhok * kernel))
    ##ex = real(fft.irfft(-1j*rhok))
    phi = N2 * np.real(np.fft.irfft(rhok * kernel))
    ##phi = real(fft.irfft(rhok))
    # phi = real(fft.irfft(rhok*kernel))
    return ex, phi


def advect_x(f, meshx, meshv, deltat):
    """Advection in x"""
    for j in range(meshv.ncells):
        alpha = meshv.getx(j) * deltat
        f[j, :] = interpSpline(3, f[j, :], meshx, alpha)
    return f


def advect_v(f, ex, meshx, meshv, deltat):
    """advection in the v direction"""
    for i in range(meshx.ncells):
        alpha = ex[i] * deltat
        f[:, i] = interpSpline(3, f[:, i], meshv, alpha)
    return f


meshv = UniformMesh(-8., 8., 64)
meshx = UniformMesh((-1) * np.pi, np.pi, 64)
X, V = np.meshgrid(meshx.getpoints(), meshv.getpoints())
eps = 0.1
# k = 2*pi/(meshx.xmax-meshx.xmin)
b = 1 / T
m = mag(b)
# m=0.0
w = np.sqrt(m)
eta = np.exp(-b * (((V ** 2) / 2) - m * np.cos(X)))
a = masse / Norm(eta, meshx, meshv)
f0 = a * np.exp(-b * (((V ** 2) / 2) - m * np.cos(X))) * (1 + eps * np.cos(X))
f = f0.copy()
ex, phi = poisson(f, meshx, meshv)

# time loop
nbiter = 10000
deltat = 0.1
T = np.zeros(nbiter)

start = time.time()

for n in range(nbiter):
    # g = abs(fft.fftshift(fft.fft2(f)))
    ##pl.clf()
    # pl.axis([2*pi - 1.0, 2*pi + 1.0, -8.0, 8.0])
    ##pl.title('Vlasov-HMF, iter='+str(n+1))
    # pl.xlabel('Velocity V')
    # pl.ylabel('Space variable X')
    # fig = pl.figure()
    # fig.add_subplot(111, projection='3d')
    # ax = Axes3D(fig)
    # ax = fig.gca(projection='3d')
    # ax.plot_surface(X,V,f)
    # lot(X,V,g)
    ##fmin=f.min()
    ##fmax=f.max()
    ##niveaux = fmin + 0.00001 + np.arange(100)*(fmax - fmin)/100.0
    ##pl.contour(X,V,f,niveaux)
    ##pl.grid()
    # extent=[0, 1, 0, 1])
    # pl.savefig('fig_'+str(n)+'.png')
    ##pl.show()
    ##pause(0.01)

    """Computation of the desired gamma for f at time n*deltat"""
    test1 = f * np.cos(X)
    # test2=f*sin(X)
    gamma1 = Norm(test1, meshx, meshv)
    # gamma2=Norm(test2,meshx,meshv)
    T[n] = gamma1

    """Vlasov-HMF Strang"""
    f = advect_x(f, meshx, meshv, deltat / 2)
    ex, phi = poisson(f, meshx, meshv)
    f = advect_v(f, ex, meshx, meshv, deltat)
    f = advect_x(f, meshx, meshv, deltat / 2)

    # free transport
    # f=advect_x(f,meshx,meshv,deltat)

    # begin Fourier cut-off
    ##filt = fft.rfft(f[:,0])
    ##nv = size(filt)
    ##nc = 3*nv/4
    ##nd = nv
    ##cutoff = zeros(nv - nc,dtype = float)
    ##inter = (arange(nc,nd) -  float(nc))/(float(nd) - float(nc))
    # cutoff[0:nd -nc] = 2*inter**3 - 3*inter**2 + 1
    # cutoff[0:nd -nc] = 6*inter**2 - 6*inter
    # for i in range(meshx.ncells):
    # filt = fft.rfft(f[:,i])
    # filt[(nc):nv] = (filt[(nc):nv])*exp(deltat*cutoff)
    # f[:,i] = fft.irfft(filt)
    # end Fourier cut-off

end = time.time()

print("Elapsed time ", end - start)

"""Substracting from gamma its long time average """
Final1 = f * np.cos(X)
Gamma1 = Norm(Final1, meshx, meshv)
for n in range(nbiter):
    T[n] = T[n] - Gamma1

mesht = UniformMesh(deltat, (nbiter - 1) * deltat, nbiter)
t = mesht.getpoints()

"""Power law fitting"""
nd = 1001  # indice pour choix de la constante
td = (nd - 1) * deltat  # temps correspondant
guess = np.array([abs(T[nd]) * (td ** 3), 3.0])
param, variance = curve_fit(power, t, T, guess, method='lm')
Fit = param[0] / (t ** param[1])
print(param)
print(variance)

figure1 = plt.figure()
plt.plot(t, np.abs(T), 'b', label="|C[f](t)-C[f][T]|")
plt.plot(t, Fit, 'r', linewidth=2, label="Envelop Fitting")
plt.xscale('symlog')
plt.yscale('symlog', linthreshy=10 ** (-12))
plt.grid(True)
plt.xlabel("t")
plt.ylabel("|C[f](t)-C[f][T]|, T=nbiter*deltat")
plt.legend()
plt.title("Decay of Gamma1 for a cosine perturbation-Envelop Fitting")
plt.savefig('fig.png')



# figure1=plt.figure()
# plt.plot(t,abs(T),label="|C[f](t)-C[f][T]|")
# plt.xscale('symlog')
# plt.yscale('symlog',linthreshy=10**(-12))
# plt.grid(True)
# plt.xlabel("t")
# plt.ylabel("|C[f](t)-C[f][T]|, T=nbiter*deltat")
# plt.legend()
# plt.title("Decay of Gamma1 for a cosine perturbation")
# plt.savefig('fig.png')
