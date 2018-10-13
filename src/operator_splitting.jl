
export @Lie, @Strang, @TripleJump, @Order6

"""

    @Lie( push_t, push_v )

    Apply the first order Lie splitting

    push_t and push_v are two function calls with
    `dt` as argument.

"""
macro Lie(push_t, push_v)
    return esc( quote
        $push_t
	$push_v
    end)
end

"""

    @Strang( push_t, push_v )

    Apply the second order Strang splitting

    push_t and push_v are two function calls with
    `dt` as argument.

"""
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

"""

    @TripleJump( push_t, push_v )

    Apply the fourth order Triple Jump splitting

    push_t and push_v are two function calls with
    `dt` as argument.

"""
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

"""

    @Order6( push_t, push_v )

    Apply the sixth order splitting

    push_t and push_v are two function calls with
    `dt` as argument.

"""
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
