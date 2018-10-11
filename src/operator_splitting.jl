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

export @Lie, @Strang, @TripleJump, @Order6

"""
    Apply the composition method for given number of times steps.
    	split	: OperatorSplitting object
    	dt	: time step
    	number_time_steps	: number of time steps to be performed
"""
macro Lie(push_t, push_v)
    return esc( quote
        $push_t
	$push_v
    end)
end

macro Strang(push_t, push_v)
    return esc(quote    
        local full_dt = dt
	dt = 0.5full_dt
        $push_t
	dt = full_dt
        $push_v
	dt = 0.5full_dt
        $push_t
	dt = full_dt
    end)
end

macro TripleJump(push_t, push_v)
    return esc(quote    
        local full_dt = dt
        dt =  0.675603595979829full_dt
        $push_t
        dt =  1.351207191959658full_dt
        $push_v
        dt = -0.17560359597982855full_dt
        $push_t
        dt = -1.702414383919315full_dt
        $push_v
        dt = -0.17560359597982855full_dt
        $push_t
        dt =  1.351207191959658full_dt
        $push_v
        dt =  0.675603595979829full_dt
        $push_t
	dt = full_dt
    end)
end

macro Order6(push_t, push_v)
    return esc(quote    
        local full_dt = dt
        dt =  0.0414649985182624full_dt
        $push_t
        dt =  0.123229775946271full_dt
        $push_v
        dt =  0.198128671918067full_dt
        $push_t
        dt =  0.290553797799558full_dt
        $push_v
        dt = -0.0400061921041533full_dt
        $push_t
        dt = -0.127049212625417full_dt
        $push_v
        dt =  0.0752539843015807full_dt          
        $push_t
        dt = -0.246331761062075full_dt
        $push_v
        dt = -0.0115113874206879full_dt
        $push_t
        dt =  0.357208872795928full_dt
        $push_v
        dt =  0.23666992478693111full_dt
        $push_t
        dt =  0.20477705429147008full_dt
        $push_v
        dt =  0.23666992478693111full_dt
        $push_t
        dt =  0.357208872795928full_dt
        $push_v
        dt = -0.0115113874206879full_dt
        $push_t
        dt = -0.246331761062075full_dt
        $push_v
        dt =  0.0752539843015807full_dt          
        $push_t
        dt = -0.127049212625417full_dt
        $push_v
        dt = -0.0400061921041533full_dt
        $push_t
        dt =  0.290553797799558full_dt
        $push_v
        dt =  0.198128671918067full_dt
        $push_t
        dt =  0.123229775946271full_dt
        $push_v
        dt =  0.0414649985182624full_dt
        $push_t
	dt = full_dt
    end)
end
