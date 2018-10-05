export PeriodicDomain


"""

    RectMesh1D1V( xmin, xmax, nx, vmin, vmax, nv)
    RectMesh1D1V( x, v)

    Regular cartesian 2D mesh for 1D1V simulation

"""
struct PeriodicDomain
    
    left   :: Float64
    right  :: Float64
    ncells :: Int
    points :: Vector{Float64}
    delta  :: Float64
    
    function PeriodicDomain(points::Vector{Float64})

        left, right, ncells = points[1], points[end], size(x)-1
	  delta = (right - left) / ncells
        new(left, right, ncells, points, delta)

    end

    function PeriodicDomain(left, right, ncells::Int)

        points = range(left, stop=right, length=ncells)[1:end-1]
	  delta  = (right - left) / ncells
        new(left, right, ncells, points, delta)

    end

end


import Base.*

Base.:*(x::PeriodicDomain, y::PeriodicDomain) = begin

    RectMesh1D1V( x.points, y.points )

end
