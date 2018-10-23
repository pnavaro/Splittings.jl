
export @Lie, @Strang, @TripleJump, @Order6, @Magic, @SuperMagic

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

    @Magic( push_t, push_v )

    Apply the second order modified Strang splitting

    with tan(dt/2) instead of dt/2 for push_t and sin(dt)
    instead of dt for push_v
    
    push_t and push_v are two function calls with
    `dt` as argument.

"""
macro Magic(push_t, push_v)
    return esc(quote    
        local full_dt = dt
	  dt = tan(0.5full_dt)
          $push_t
	  dt = sin(full_dt)
          $push_v
	  dt = tan(0.5full_dt)
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

macro SuperMagic(push_t, push_v)
    return esc(quote

        local ct = cos(dt)
	local st = sin(dt)
        local α  = (1 - ct)/st
	local β  = (dt - st*ct)/(24*(ct^2-1))
	local γ  = ((-2dt + 3st) * ct + dt - 2st)/(6*(ct-1)*(ct+1)^2)
	local μ  = st / 2 
	local ν  = - (2dt + 3st) * ct^2 / 192 - (3dt + 5st)*ct/192 - dt/192 + 7st/96
        local full_dt = dt

	dt = 1/full_dt
	v  .= α .* v .+ β * v.^3
	$push_t
	dt = 1/full_dt
	e  .= μ * e .+ ν * x.^3
	$push_v
	dt = 1/full_dt
	v  .= γ * v.^3
	$push_t
	dt = 1/full_dt
	e  .= μ * e .+ ν * x.^3
	$push_v
	dt = 1/full_dt
	v  .= α .* v .+ β * v.^3
	$push_t


    end)
end
