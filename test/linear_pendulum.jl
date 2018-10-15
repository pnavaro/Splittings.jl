import Splittings: @Lie, @Strang, @TripleJump, @Order6
"""

    Implements split operators for linear pendulum 

    Solve linear pendulum problem: ``\\frac{delta1}{dt} = v ``, 

    ``\\frac{delta2}{dt} = - ω^2 x``. The exact solution is
    ``x(t)=  x(0)   cos(ωt) + \frac{v(0)}{ω}sin(ωt),```
    ``v(t)= -x(0) ω sin(ωt) + v(0)cos(ωt) ``

"""
function exact( t_final, ω, x0, v0 )
    x = x0 * cos( ω * t_final ) + (v0 / ω) * sin( ω * t_final )
    v = -x0 * ω * sin( ω* t_final ) + v0 * cos( ω * t_final )
    return x, v
end

function check_order(do_split_steps::Function, 
		     steps_fine::Int64, 
		     expected_order::Int )

    ω  = 2.0
    x0 = 1.0       # initial x  for order checking
    v0 = 2.0       # initial v  for order checking 
    t_final = 1.0  # final time for order checking 

    x_exact, v_exact = exact( t_final, ω, x0, v0 )

    # do iterations with smallest time step
    dt = t_final/steps_fine
    number_time_steps = steps_fine

    x, v = do_split_steps((x0, v0), dt, number_time_steps, ω)
  
    # compute  mean square error
    error0 = sqrt( (x - x_exact)^2 + (v - v_exact)^2 )

    # do iterations with middle time step
    dt = 2dt
    number_time_steps = number_time_steps ÷ 2

    x, v = do_split_steps((x0, v0), dt, number_time_steps, ω)
 
    # compute mean square error
    error1 = sqrt( (x - x_exact)^2 + (v - v_exact)^2 )
  
    # do iterations with largest time step
    dt = 2dt
    number_time_steps = number_time_steps ÷ 2

    x, v = do_split_steps((x0, v0), dt, number_time_steps, ω)
  
    # compute  mean square error
    error2 = sqrt( (x - x_exact)^2 + (v - v_exact)^2 )

    # compute order
    order1 = log(error1/error0)/log(2.0)
    order2 = log(error2/error1)/log(2.0)

    if (    (abs(order1-expected_order) > 5.e-2) 
         || (abs(order2-expected_order) > 5.e-2))

       println( "error coarse = $error2")
       println( "error middle = $error1")
       println( "      order (coarse/middle) = $order2")
       println( "error fine   = $error0")
       println( "      order (middle/fine)   = $order1")

       return false

    end

    true

end 

" First operator of splitting for linear pendulum"
function push_x!( x, v, dt ) 
    x .= x .+ v * dt 
end

" Second operator of splitting for linear pendulum"
function push_v!( x, v, dt, ω) 
    v .= v .- ω^2 * x * dt 
end
  
@testset "Splitting Operators macros" begin

    function do_steps_lie_tv( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
	x = zeros(Float64, 1)
	v = zeros(Float64, 1)
	x[1], v[1] = start
        for i = 1:nsteps
            @Lie push_x!(x,v,dt) push_v!(x,v,dt,ω)
        end
	x[1], v[1]
    end
    
    @test check_order(do_steps_lie_tv, 200, 1)
    
    function do_steps_lie_vt( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
	x = zeros(Float64, 1)
	v = zeros(Float64, 1)
	x[1], v[1] = start
        for i = 1:nsteps
	    @Lie push_v!(x,v,dt,ω) push_x!(x,v,dt)
        end
	x[1], v[1]
    end
    
    @test check_order(do_steps_lie_vt, 200, 1)
    
    function do_steps_strang_tvt( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
	x = zeros(Float64, 1)
	v = zeros(Float64, 1)
	x[1], v[1] = start
        for i = 1:nsteps
            @Strang push_x!(x,v,dt) push_v!(x,v,dt,ω)
        end
	x[1], v[1]
    end
    
    @test check_order(do_steps_strang_tvt, 100, 2)
    
    function do_steps_strang_vtv( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
	x = zeros(Float64, 1)
	v = zeros(Float64, 1)
	x[1], v[1] = start
        for i = 1:nsteps
            @Strang push_v!(x,v,dt,ω) push_x!(x,v,dt) 
        end
	x[1], v[1]
    end
    
    @test check_order(do_steps_strang_vtv, 100, 2)

    function do_steps_triple_jump_tvt( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
	x = zeros(Float64, 1)
	v = zeros(Float64, 1)
	x[1], v[1] = start
        for i = 1:nsteps
            @TripleJump push_x!(x,v,dt) push_v!(x,v,dt,ω)
        end
	x[1], v[1]
    end
    
    @test check_order(do_steps_triple_jump_tvt, 64, 4)

    function do_steps_triple_jump_vtv( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
	x = zeros(Float64, 1)
	v = zeros(Float64, 1)
	x[1], v[1] = start
        for i = 1:nsteps
            @TripleJump push_v!(x,v,dt,ω) push_x!(x,v,dt) 
        end
	x[1], v[1]
    end
    
    @test check_order(do_steps_triple_jump_vtv, 64, 4)

    function do_steps_order6_tvt( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
	x = zeros(Float64, 1)
	v = zeros(Float64, 1)
	x[1], v[1] = start
        for i = 1:nsteps
            @Order6 push_x!(x,v,dt) push_v!(x,v,dt,ω)
        end
	x[1], v[1]
    end
    
    @test check_order(do_steps_order6_tvt, 20, 6)

    function do_steps_order6_vtv( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
	x = zeros(Float64, 1)
	v = zeros(Float64, 1)
	x[1], v[1] = start
        for i = 1:nsteps
            @Order6 push_v!(x,v,dt,ω) push_x!(x,v,dt) 
        end
	x[1], v[1]
    end
    
    @test check_order(do_steps_order6_vtv, 20, 6)

end
