using IntervalSets

export PeriodicDomain

"""

    PeriodicDomain( left, right, ncells)
    PeriodicDomain( left..right, ncells)
    PeriodicDomain( points )

    Regular cartesian 2D mesh for 1D1V simulation

"""
mutable struct PeriodicDomain
    
    left   :: Float64
    right  :: Float64
    ncells :: Int
    points :: Vector{Float64}
    delta  :: Float64
    bc     :: Symbol
    axis   :: Int
    
    function PeriodicDomain(points::Vector{Float64})

        left, right, ncells = points[1], points[end], size(x)-1
	  delta = (right - left) / ncells
        new(left, right, ncells, points, delta)
	  bc = :periodic

    end

    function PeriodicDomain(left, right, ncells::Int)

        points = range(left, stop=right, length=ncells+1)[1:end-1]
	  delta  = (right - left) / ncells
	  new(left, right, ncells, points, delta, :periodic, 1)

    end

    function PeriodicDomain(interval::ClosedInterval{Int64}, ncells::Int)

	  left  = interval.left
	  right = interval.right
        points = range(left, stop=right, length=ncells+1)[1:end-1]
	  delta  = (right - left) / ncells
        new(left, right, ncells, points, delta, :periodic, 1)

    end

    function PeriodicDomain(s::StepRangeLen)

	  left   = s.offset
	  ncells = s.len-1
	  delta  = s.step
	  right  = left + ncells * delta 
	  points = collect(s)
        new(left, right, ncells, points, delta, :periodic, 1)

    end

end


import Base.*

Base.:*(x::PeriodicDomain, y::PeriodicDomain) = begin

    if ( x.axis >= y.axis )
        RectMesh1D1V( y.points, x.points )
    else
        RectMesh1D1V( x.points, y.points )
    end

end

import LinearAlgebra

import LinearAlgebra.transpose


LinearAlgebra.:transpose(d::PeriodicDomain) = begin
    
    d.axis = (d.axis+2)%2+1
    d

end
