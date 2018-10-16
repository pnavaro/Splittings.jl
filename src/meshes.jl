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
