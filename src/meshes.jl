"""

    RectMesh1D1V( x1min, x1max, n1, x2min, x2max, n2)
    RectMesh1D1V( x, v)

    Regular cartesian 2D mesh for 1D1V simulation

"""
struct RectMesh1D1V
    
    x1min :: Float64
    x1max :: Float64
    n1   :: Int
    x2min :: Float64
    x2max :: Float64
    n2   :: Int
    x    :: Vector{Float64}
    v    :: Vector{Float64}
    delta1   :: Float64
    delta2   :: Float64
    bc   :: Dict{Symbol, Symbol}
    
    function RectMesh1D1V(x::Vector{Float64}, v::Vector{Float64})

        n1 = size(x)[1]
        n2 = size(v)[1]
        x1min, x1max = x[1], x[end]
        x2min, x2max = v[1], v[end]
        delta1 = (x1max - x1min) / (n1-1)
        delta2 = (x2max - x2min) / (n2-1)
        bc = Dict(:x=>:periodic,:v=>:periodic)
        new(x1min, x1max, n1, x2min, x2max, n2, x, v, delta1, delta2, bc)

    end

    function RectMesh1D1V(x1min, x1max, n1::Int,  
                          x2min, x2max, n2::Int)

        x  = range(x1min, stop=x1max, length=n1)
        v  = range(x2min, stop=x2max, length=n2)
        delta1 = (x1max - x1min) / (n1-1)
        delta2 = (x2max - x2min) / (n2-1)
        bc = Dict(:x=>:periodic,:v=>:periodic)
        new(x1min, x1max, n1, x2min, x2max, n2, x, v, delta1, delta2, bc)

    end

    function RectMesh1D1V(x::StepRangeLen, v::StepRangeLen)

	  x1min = x.offset
	  n1   = x.len
	  delta1   = x.step
	  x1max = x1min + (n1-1) * delta1 

	  x2min = v.offset
	  n2   = v.len
	  delta2   = v.step
	  x2max = x2min + (n2-1) * delta2 

	  bc = Dict(:x=>:periodic,:v=>:periodic)
        new(x1min, x1max, n1, x2min, x2max, n2, x, v, delta1, delta2, bc)

    end

end


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

       if (endpoint)
           points = range(start, stop=stop, length=length)
       else
           points = range(start, stop=stop, length=length+1)[1:end-1]
       end
       step = points[2]-points[1]

       new( start, stop, length, step, points, endpoint)

   end

end
