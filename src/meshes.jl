"""

    RectMesh1D1V( xmin, xmax, nx, vmin, vmax, nv)
    RectMesh1D1V( x, v)

    Regular cartesian 2D mesh for 1D1V simulation

"""
struct RectMesh1D1V
    
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
    bc   :: Dict{Symbol, Symbol}
    
    function RectMesh1D1V(x::Vector{Float64}, v::Vector{Float64})

        nx = size(x)[1]
        nv = size(v)[1]
        xmin, xmax = x[1], x[end]
        vmin, vmax = v[1], v[end]
	dx = (xmax - xmin) / (nx-1)
	dv = (vmax - vmin) / (nv-1)
	bc = Dict(:x=>:periodic,:v=>:periodic)
        new(xmin, xmax, nx, vmin, vmax, nv, x, v, dx, dv, bc)

    end

    function RectMesh1D1V(xmin, xmax, nx::Int,  
                          vmin, vmax, nv::Int)

        x  = range(xmin, stop=xmax, length=nx)
        v  = range(vmin, stop=vmax, length=nv)
	dx = (xmax - xmin) / (nx-1)
	dv = (vmax - vmin) / (nv-1)
	bc = Dict(:x=>:periodic,:v=>:periodic)
        new(xmin, xmax, nx, vmin, vmax, nv, x, v, dx, dv, bc)

    end

    function RectMesh1D1V(x::StepRangeLen, v::StepRangeLen)

	xmin = x.offset
	nx   = x.len
	dx   = x.step
	xmax = xmin + (nx-1) * dx 

	vmin = v.offset
	nv   = v.len
	dv   = v.step
	vmax = vmin + (nv-1) * dv 

	bc = Dict(:x=>:periodic,:v=>:periodic)
        new(xmin, xmax, nx, vmin, vmax, nv, x, v, dx, dv, bc)

    end

end

"""

    RectMesh1D(xmin, xmax, nx)

    1D uniform mesh data

"""
struct RectMesh1D

   xmin     :: Float64
   xmax     :: Float64
   nx       :: Int
   dx       :: Float64
   x        :: Vector{Float64}
   endpoint :: Bool

   function RectMesh1D(xmin, xmax, nx::Int; endpoint=true)

       dx = (xmax - xmin) / nx
       if (endpoint)
           x  = range(xmin, stop=xmax, length=nx)
       else
           x  = range(xmin, stop=xmax, length=nx+1)[1:end-1]
       end

       new( xmin, xmax, nx, dx, x, endpoint)

   end

end
