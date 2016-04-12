### substitute ########################################

substitute(ex::TimeVariable, this::TimeVariable, by::TimeExpression) = (ex==this ? by : ex)

function substitute(ex::TimeLinearCombination, this::TimeVariable, by::TimeExpression)
    TimeLinearCombination(Tuple{TimeExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

substitute(ex::SpaceVariable, this::SpaceVariable, by::SpaceExpression) = (ex==this ? by : ex)
substitute(ex::SpaceVariable, this::TimeVariable, by::TimeExpression) = ex 

function substitute(ex::SpaceLinearCombination, this::SpaceVariable, by::SpaceExpression)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

function substitute(ex::SpaceLinearCombination, this::TimeVariable, by::TimeExpression)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

function substitute(ex::AutonomousFunctionExpression, this::SpaceVariable, by::SpaceExpression)
    ex.fun(substitute(ex.x, this, by), [substitute(ex.c, this, by) for x in ex.d_args]...)
end

function substitute(ex::FlowExpression, this::SpaceVariable, by::SpaceExpression)
    E(ex.fun, ex.t, substitute(ex.x, this, by), [substitute(ex.c, this, by) for x in ex.d_args]...)
end

# TODO: to be completed

### differential #######################################

differential(ex::SpaceVariable, with_respect_to::SpaceVariable, applied_to::SpaceExpression) =
      (ex==with_respect_to ? applied_to : x_zero)

function differential(ex::SpaceLinearCombination, with_respect_to::SpaceVariable, applied_to::SpaceExpression)   
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(differential(x, with_respect_to, applied_to), c) 
        for (x, c) in ex.terms])
end

function differential(ex::AutonomousFunctionExpression, with_respect_to::SpaceVariable, 
                      applied_to::SpaceExpression)    
    terms = Array{SpaceExpression}(length(ex.d_args)+1)
    dx = differential(ex.x, with_respect_to, applied_to)   
    terms[1] = AutonomousFunctionExpression(ex.fun, ex.x, vcat(ex.d_args, dx)...)
    for i=1:length(ex.d_args)
          d_args = copy(ex.d_args)
          d_args[i] = differential(ex.d_args[i], with_respect_to, applied_to)
          terms[i+1] = AutonomousFunctionExpression(ex.fun, ex.x, d_args...)
    end      
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(x, 1) for x in terms])
end

function differential(ex::FlowExpression, with_respect_to::SpaceVariable, applied_to::SpaceExpression)    
    terms = Array{SpaceExpression}(length(ex.d_args)+1)
    dx = differential(ex.x, with_respect_to, applied_to)   
    terms[1] = FlowExpression(ex.fun, ex.t, ex.x, ex.dt_order, vcat(ex.d_args, dx)...)
    for i=1:length(ex.d_args)
          d_args = copy(ex.d_args)
          d_args[i] = differential(ex.d_args[i], with_respect_to, applied_to)
          terms[i+1] = FlowExpression(ex.fun, ex.t, ex.x, ex.dt_order, d_args...)
    end      
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(x, 1) for x in terms])
end


### t_derivative #######################################

t_derivative(x::SpaceVariable, with_respect_to::TimeVariable; flag::Bool=false) = x_zero 
    # time derivative of a space variable...

function t_derivative(ex::SpaceLinearCombination, with_respect_to::TimeVariable; flag::Bool=false)   
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(t_derivative(x, with_respect_to, flag=flag), c) 
        for (x, c) in ex.terms])
end

function t_derivative(ex::AutonomousFunctionExpression, with_respect_to::TimeVariable; flag::Bool=false)    
    terms = Array{SpaceExpression}(length(ex.d_args)+1)
    dx = t_derivative(ex.x, with_respect_to, flag=flag)    
    terms[1] = AutonomousFunctionExpression(ex.fun, ex.x, vcat(ex.d_args, dx)...)
    for i=1:length(ex.d_args)
          d_args = copy(ex.d_args)
          d_args[i] = t_derivative(ex.d_args[i], with_respect_to, flag=flag)
          terms[i+1] = AutonomousFunctionExpression(ex.fun, ex.x, d_args...)
    end      
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(x, 1) for x in terms])
end

global x_var = SpaceVariable("x_var")

function t_derivative(ex::FlowExpression, with_respect_to::TimeVariable; flag::Bool=false)    
    f = coefficient(ex.t, with_respect_to)
    dx = t_derivative(ex.x, with_respect_to, flag=flag)    
    terms = Array{SpaceExpression}(length(ex.d_args)+2)
    if f != 0
        if !flag
            t = AutonomousFunctionExpression(ex.fun, FlowExpression(ex.fun, ex.t, x_var, 0))
        else
            t = FlowExpression(ex.fun, ex.t, x_var, 0, AutonomusFunctionExpression(ex.fun, x_var))
        end
        for d_arg in ex.d_args
            t = differential(t, x_var, d_arg)
        end
        terms[1] = f*substitute(t, x_var, ex.x) 
    else
        terms[1] = x_zero
    end
    terms[2] = FlowExpression(ex.fun, ex.t, ex.x, ex.dt_order, vcat(ex.d_args, dx)...)
    for i=1:length(ex.d_args)
          d_args = copy(ex.d_args)
          d_args[i] = t_derivative(ex.d_args[i], with_respect_to, flag=flag)
          terms[i+2] = FlowExpression(ex.fun, ex.t, ex.x, ex.dt_order, d_args...)
    end      
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(x, 1) for x in terms])
end

