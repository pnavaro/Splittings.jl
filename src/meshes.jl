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

export UniformMesh

"""

    UniformMesh(start, stop, length)

    1D uniform mesh data.

    length   : number of points
    length-1 : number of cells

    To remove the last point, set endpoint=false

"""
struct UniformMesh

   start    :: Float64
   stop     :: Float64
   length   :: Int
   step     :: Float64
   points   :: Vector{Float64}
   endpoint :: Bool

   function UniformMesh(start, stop, length::Int; endpoint=true)

       step = (stop - start) / (length-1)
       if (endpoint)
           points = range(start, stop=stop, length=length)
       else
           points = range(start, stop=stop, length=length+1)[1:end-1]
       end

       new( start, stop, length, step, points, endpoint)

   end

   function UniformMesh(xrange::StepRangeLen)

       step     = xrange.step
       start    = xrange.offset
       length   = xrange.len
       stop     = start + (length-1) * step
       points   = collect(xrange)
       endpoint = true

       new( start, stop, length, step, points, endpoint)

   end

end
