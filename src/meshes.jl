"""

   Simple regular cartesian 2D mesh

"""
struct Mesh
    nx   :: Int
    ny   :: Int
    xmin :: Float64
    xmax :: Float64
    ymin :: Float64
    ymax :: Float64
end

"""

   Regular cartesian 2D mesh for 1D1V simulation

"""
struct Mesh1D1V
    
    xmin :: Float64
    xmax :: Float64
    nx   :: Int
    vmin :: Float64
    vmax :: Float64
    nv   :: Int
    x    :: Vector{Float64}
    v    :: Vector{Float64}
    dx   :: Float64
    dv   :: Float64
    
    function Mesh1D1V(xmin, xmax, nx, vmin, vmax, nv)
        x = range(xmin, stop=xmax, length=nx)
        v = range(vmin, stop=vmax, length=nv)
	dx = (xmax - xmin) / (nx-1)
	dv = (vmax - vmin) / (nv-1)
        new(xmin, xmax, nx, vmin, vmax, nv, x, v, dx, dv)
    end
end

"""

   2D rectangular cartesian mesh parameters

"""
struct RectMesh2D

    xmin :: Float64
    xmax :: Float64
    nx   :: Int
    dx   :: Float64
    ymin :: Float64
    ymax :: Float64
    ny   :: Int
    dy   :: Float64

    function RectMesh2D(xmin, xmax, nx, ymin, ymax, ny)
       dx = (xmax - xmin) / nx
       dy = (ymax - ymin) / ny
       new( xmin, xmax, nx, dx, ymin, ymax, ny, dy)
    end

end

"""

   1D uniform mesh data

"""
struct UniformMesh
   xmin  :: Float64
   xmax  :: Float64
   nx    :: Int
   dx    :: Float64
   x     :: Vector{Float64}
   function UniformMesh(xmin, xmax, nx)
      dx = (xmax - xmin) / nx
      x  = range(xmin, stop=xmax, length=nx+1)[1:end-1]
      new( xmin, xmax, nx, dx, x)
   end
end
