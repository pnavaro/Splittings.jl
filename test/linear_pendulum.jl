"""

    Implements split operators for linear pendulum 

    Solve linear pendulum problem: ``\\frac{dx}{dt} = v ``, 

    ``\\frac{dv}{dt} = - ω^2 x``. The exact solution is
    ``x(t)=  x(0)   cos(ωt) + \frac{v(0)}{ω}sin(ωt),```
    ``v(t)= -x(0) ω sin(ωt) + v(0)cos(ωt) ``

"""
abstract type LinearPendulum end


" Exact solution of linear pendulum "
function exact( t_final, ω, x0, v0 )
    x = x0 * cos( ω * t_final ) + (v0 / ω) * sin( ω * t_final )
    v = -x0 * ω * sin( ω* t_final ) + v0 * cos( ω * t_final )
    return x, v
end

function check_order(do_split_steps::Function, 
		     steps_fine::Int64, 
		     expected_order::Int )

    ω  = 2.0
    x0 = 1.0       # initial x for order checking
    v0 = 2.0       # initial v for order checking 
    t_final = 1.0  # final time for order checking 

    x_exact, v_exact = exact( t_final, ω, x0, v0 )

    # do iterations with smallest time step
    dt = t_final/steps_fine
    number_time_steps = steps_fine

    x, v = do_split_steps((x0, v0), dt, number_time_steps, ω)
  
    # compute  mean square error
    error0 = sqrt( (x - x_exact)^2 + (v - v_exact)^2 )

    # do iterations with middle time step
    dt = 2*dt
    number_time_steps = number_time_steps ÷ 2

    x, v = do_split_steps((x0, v0), dt, number_time_steps, ω)
 
    # compute mean square error
    error1 = sqrt( (x - x_exact)^2 + (v - v_exact)^2 )
  
    # do iterations with largest time step
    dt = 2*dt
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
push_x( x, v, dt ) = x + v * dt

" Second operator of splitting for linear pendulum"
push_v( x, v, dt, ω) = v - ω^2 * x * dt
  
@testset "Splitting Operators macros" begin

    function do_steps( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
        x, v = start
        for i = 1:nsteps
            @LieTV push_x(x,v,dt) push_v(x,v,dt,ω)
        end
        x, v
    end
    
    @test check_order(do_steps, 200, 1)
    
    function do_steps( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
        x, v = start
        for i = 1:nsteps
            @LieVT push_x(x,v,dt) push_v(x,v,dt,ω)
        end
        x, v
    end
    
    @test check_order(do_steps, 200, 1)
    
    function do_steps( start::Tuple{Float64,Float64},
    		   dt::Float64, nsteps::Int64, ω::Float64 )
        x, v = start
        for i = 1:nsteps
            @StrangTVT push_x(x,v,dt) push_v(x,v,dt,ω)
        end
        x, v
    
    end
    
    @test check_order(do_steps, 100, 2)
    
    #function do_steps( start::Tuple{Float64,Float64},
    #		   dt::Float64, nsteps::Int64, ω::Float64 )
    #    x, v = start
    #    for i = 1:nsteps
    #        @StrangVTV push_x(x,v,dt) push_v(x,v,dt,ω)
    #    end
    #    x, v
    #end
    #
    #@test check_order(do_steps, 100, 2)

end

#tests = [(lie_tv,200, 1), (lie_vt,200, 1), (strang_tvt,100, 2),
#         (strang_vtv,100, 1), (triple_jump_tvt, 64, 4), 
#         (triple_jump_vtv,64,4), (order6_tvt,20,4),
#         (order6_vtv,20, 6)] 
#
#for (method, nsteps, expected_order) in test
#
#  if (check_order(method, nsteps, expected_order))
#     println("Test with $method passed ")
#  else
#     println("Problem with order of $method ")
#  end if
#
#
#end 

