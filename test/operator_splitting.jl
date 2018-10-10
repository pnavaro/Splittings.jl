"""
The operator splitting module provides a generic implementation of composition algorithms 
of different order for two operators T and V for solving `` \\frac{dU}{dt} = (T+V) U ``. 
The solution on one time step can be written `` U(Δt) = \\mathcal{S}_{T+V} U(0) ``. 
The composition algorithm consists in successive solutions of the split equations 
`` \\frac{dU}{dt} = T U `` and `` \\frac{dU}{dt} = V U ``. 
Alternating the two reduced solution operators `` \\mathcal{S}_{T}  ``
and `` \\mathcal{S}_{V} `` with adequately chosen time increments yields arbitrary 
order in time for the full solution.

The application of an operator splitting method to a concrete problem is done
by extending the sll_t_operator_splitting splitting base class by a new type
containing on the one hand the data on which the operators act
and a specific implementation of the two operators
"""
abstract type OperatorSplitting end

struct LieTV <: OperatorSplitting

    steps   :: Vector{Float64}
    push_t! :: Function
    push_v! :: Function
    p       :: Any

    function LieTV( push_t!, push_v! ; p=0.0)
        steps = Float64[]
        push!(steps, 1.)
        push!(steps, 1.)
        new( steps, push_t!, push_v!, p )
    end

end

struct LieVT <: OperatorSplitting

    steps   :: Vector{Float64}
    push_t! :: Function
    push_v! :: Function
    p       :: Any

    function LieTV( push_t!, push_v! ; p= 0.0)
        steps = Float64[]
        push!(steps, 1.)
        push!(steps, 1.)
        new( steps, push_t!, push_v!, p )
    end

end


"""
    Apply the composition method for given number of times steps.
    	split	: OperatorSplitting object
    	dt	: time step
    	number_time_steps	: number of time steps to be performed
"""
function do_split_steps(split::LieTV, start::Tuple{Float64, Float64}, 
    dt::Float64, number_time_steps)

    x, v = start
    for i = 1:number_time_steps
        split.push_t!(x, v, dt, split.p)
        split.push_v!(x, v, dt, split.p)
    end
    x, v

end

function do_split_steps(split::LieVT, start::Tuple{Float64, Float64},
    dt::Float64, number_time_steps)

    x, v = start
    for i = 1:number_time_steps
        split.push_v!(x, v, dt, split.p)
        split.push_t!(x, v, dt, split.p)
    end
    x, v

end

#struct StrangTVT
#    
#       split%nb_split_step = 3
#       split%split_begin_T = .true.
#       push!(split_step, 0.5_f64
#       push!(split_step, 1._f64
#       push!(split_step, split%split_step(1)
#end
#
#struct StrangVTV
#    
#       split%nb_split_step = 3
#       split%split_begin_T = .false.
#       push!(split_step(1) = 0.5_f64
#       push!(split_step(2) = 1._f64
#       push!(split_step(3) = split%split_step(1)
#end
#
#struct TripleJumpTVT
#
#       split%nb_split_step = 7
#       split%split_begin_T = .true.
#       push!(split_step(1) = 0.675603595979829_f64
#       push!(split_step(2) = 1.351207191959658_f64
#       push!(split_step(3) = -0.17560359597982855_f64
#       push!(split_step(4) = -1.702414383919315_f64
#       push!(split_step(5) = split%split_step(3)
#       push!(split_step(6) = split%split_step(2)
#       push!(split_step(7) = split%split_step(1)
#
#end
#
#struct TripleJumpVTVT
#
#       split%nb_split_step = 7
#       split%split_begin_T = .false.
#       push!(split_step(1) = 0.675603595979829_f64
#       push!(split_step(2) = 1.351207191959658_f64
#       push!(split_step(3) = -0.17560359597982855_f64
#       push!(split_step(4) = -1.702414383919315_f64
#       push!(split_step(5) = split%split_step(3)
#       push!(split_step(6) = split%split_step(2)
#       push!(split_step(7) = split%split_step(1)
#
#end
#
#struct Order6TVT
#       split%nb_split_step = 23
#       split%split_begin_T = .true.
#       push!(split_step(1) = 0.0414649985182624_f64
#       push!(split_step(2) = 0.123229775946271_f64
#       push!(split_step(3) = 0.198128671918067_f64
#       push!(split_step(4) = 0.290553797799558_f64
#       push!(split_step(5) = -0.0400061921041533_f64
#       push!(split_step(6) = -0.127049212625417_f64
#       push!(split_step(7) = 0.0752539843015807_f64          
#       push!(split_step(8) = -0.246331761062075_f64
#       push!(split_step(9) = -0.0115113874206879_f64
#       push!(split_step(10) = 0.357208872795928_f64
#       push!(split_step(11) = 0.23666992478693111_f64
#       push!(split_step(12) = 0.20477705429147008_f64
#       push!(split_step(13) = split%split_step(11)
#       push!(split_step(14) = split%split_step(10)
#       push!(split_step(15) = split%split_step(9)          
#       push!(split_step(16) = split%split_step(8)
#       push!(split_step(17) = split%split_step(7)
#       push!(split_step(18) = split%split_step(6)
#       push!(split_step(19) = split%split_step(5)
#       push!(split_step(20) = split%split_step(4)
#       push!(split_step(21) = split%split_step(3)
#       push!(split_step(22) = split%split_step(2)
#       push!(split_step(23) = split%split_step(1)  
#end
#
#struct Order6VTV
#
#       split%nb_split_step = 23
#       split%split_begin_T = .false.
#       push!(split_step(1) = 0.0414649985182624_f64
#       push!(split_step(2) = 0.123229775946271_f64
#       push!(split_step(3) = 0.198128671918067_f64
#       push!(split_step(4) = 0.290553797799558_f64
#       push!(split_step(5) = -0.0400061921041533_f64
#       push!(split_step(6) = -0.127049212625417_f64
#       push!(split_step(7) = 0.0752539843015807_f64          
#       push!(split_step(8) = -0.246331761062075_f64
#       push!(split_step(9) = -0.0115113874206879_f64
#       push!(split_step(10) = 0.357208872795928_f64
#       push!(split_step(11) = 0.23666992478693111_f64
#       push!(split_step(12) = 0.20477705429147008_f64
#       push!(split_step(13) = split%split_step(11)
#       push!(split_step(14) = split%split_step(10)
#       push!(split_step(15) = split%split_step(9)          
#       push!(split_step(16) = split%split_step(8)
#       push!(split_step(17) = split%split_step(7)
#       push!(split_step(18) = split%split_step(6)
#       push!(split_step(19) = split%split_step(5)
#       push!(split_step(20) = split%split_step(4)
#       push!(split_step(21) = split%split_step(3)
#       push!(split_step(22) = split%split_step(2)
#       push!(split_step(23) = split%split_step(1)          
#
#end
#
#
#"""
#
#
#"""
# Combine first and last step of symmetric operators
# if more than one time step if performed
#"""
#function do_split_steps(split::OperatorSplitting, dt::Float64, number_time_steps)
#
#    if (split%split_begin_T)
#       # Apply T first
#       # First split step
#       istep = 1 
#       operatorT(split, split%split_step(istep)*dt)
#       for i = 1:number_time_steps - 1
#          # Start with second split step as first already accounted
#          # for in last time step
#          istep = 2
#          while (istep < split.nb_split_step - 2)
#             operatorV(split, split%split_step(istep)*dt)
#             istep = istep + 1
#             operatorT(split, split%split_step(istep)*dt)
#             istep = istep + 1
#          end do
#          # Last two steps with T applied on twice the step as 
#          # it is combine with the first push of the next time step
#          operatorV(split, split%split_step(istep)*dt)
#          istep = istep + 1
#          operatorT(split, 2*split%split_step(istep)*dt)
#          istep = istep + 1
#       end
#       # Last time step of the sequence. Start with second split step
#       # as first already accounted for in last time step
#       istep = 2
#       while (istep < split%nb_split_step)
#          operatorV(split, split%split_step(istep)*dt)
#          istep = istep + 1
#          operatorT(split, split%split_step(istep)*dt)
#          istep = istep + 1
#       end
#    else     
#       # Apply V first
#       # First split step
#       istep = 1 
#       operatorV(split, split%split_step(istep)*dt)
#       for i = 1:number_time_steps - 1
#          # Start with second split step as first already accounted
#          # for in last time step
#          istep = 2
#          while (istep < split%nb_split_step - 2)
#             operatorT(split, split_step(istep)*dt)
#             istep = istep + 1
#             operatorV(split, split_step(istep)*dt)
#             istep = istep + 1
#          end
#          # Last two steps with V applied on twice the step as 
#          # it is combine with the first push of the next time step
#          operatorT(split, split.split_step(istep)*dt)
#          istep = istep + 1
#          operatorV(split, 2*split%split_step(istep)*dt)
#          istep = istep + 1
#       end
#       ! Last time step of the sequence. Start with second split step
#       ! as first already accounted for in last time step
#       istep = 2
#       while (istep < split.nb_split_step) 
#          operatorT(split, split.split_step(istep)*dt)
#          istep = istep + 1
#          operatorV(split, split.split_step(istep)*dt)
#          istep = istep + 1
#       end
#    end
#
#end 
