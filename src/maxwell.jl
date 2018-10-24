hx = hx - dt ( iffty(-jky ffty(ez) 
	+ dt^2/24( ifftx(-kx^2 fftx ( iffty ( -jky ffty(ez)))))
        + iffty(jky^3 ffty(ez))))

hy = hy + dt ( ifftx(-jkx fftx(ez)) 
	+ dt^2/24 *(iffty(-ky^2 ffty(ifftx(-jkx fftx(ez)))) 
	+ ifftx(jkx^3 fftx(ez))))

ez = ez + dt ( ifftx(-jkx fftx(hy)) - iffty(-jky ffty(hx)) 
	+ dt^2/24 ( 
	      2 iffty( -ky^2 ffty(ifftx(-jkx fftx(hy))))
	      + ifftx( jkx^3 fftx(hy)) 
	      - iffty( jky^3 ffty(hx))))


ex = ex + dt ( iffty(-jky ffty(hz)
        + dt^2/24 ( ifftx(-kx^2 fftx ( iffty ( -jky ffty(hz)))))
        + iffty(jky^3 ffty(hz))))

ey = ey - dt ( ifftx(-jkx fftx(hz)) 
	+ dt^2/24 *(iffty(-ky^2 ffty(ifftx(-jkx fftx(hz)))) 
	+ ifftx(jkx^3 fftx(hz))))

hz = hz - dt (  ifftx(-jkx fftx(ey)) 
	      - iffty(-jky ffty(ex)) 
	      - dt^2/24 ( 
	      - iffty(-ky^2 ffty(ifftx(-jkx fftx(ey))))
	      + ifftx(-kx^2 fftx(iffty(-jky ffty(ex))))
	      + ifftx(jkx^3 fftx(ey)) 
              - iffty(jky^3 ffty(ex))))

omega      =   sqrt(8*pi*pi)
bz         = - cos(2*pi*x)*cos(2*pi*y)*cos(omega*t)
ex         =   cos(2*pi*x)*sin(2*pi*y)*sin(omega*t)*2*pi/ omega
ey         = - sin(2*pi*x)*cos(2*pi*y)*sin(omega*t)*2*pi/ omega
diff(ex,t) =   diff(bz,y)
diff(ey,t) = - diff(bz,x)
diff(bz,t) =   diff(ex,y) - diff(ey,x)

omega = sqrt(2*pi*pi)
ez =  cos(pi*x)*cos(pi*y)*cos(omega*t)
hx =  cos(pi*x)*sin(pi*y)*sin(omega*t) * pi / omega
hy = -sin(pi*x)*cos(pi*y)*sin(omega*t) * pi / omega

diff(hx,t) = - diff(ez,y)
diff(hy,t) =   diff(ez,x)
diff(ez,t) =   diff(hy,x) - diff(hx,y)

  d(d2hz/dy2 + d2hz/dx2)/dy
- d(d2hz/dx2 + d2hz/dy2)/dx
  d3ex/dy3 - d3dey/dx3 + d3ex/dydx2 - d3ey/dxdy2 
