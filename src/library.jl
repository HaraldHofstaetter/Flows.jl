global x_var = SpaceVariable("x_var")

### substitute ########################################

# substitute TimeVariable by TimeExpression

substitute(ex::TimeVariable, this::TimeVariable, by::TimeExpression) = (ex==this ? by : ex)

function substitute(ex::TimeLinearCombination, this::TimeVariable, by::TimeExpression)
    TimeLinearCombination(Tuple{TimeExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

substitute(ex::SpaceVariable, this::TimeVariable, by::TimeExpression) = ex 

function substitute(ex::SpaceLinearCombination, this::TimeVariable, by::TimeExpression)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

function substitute(ex::AutonomousFunctionExpression, this::TimeVariable, by::TimeExpression)
    AutonomousFunctionExpression(ex.fun, substitute(ex.x, this, by), [substitute(x, this, by) for x in ex.d_args]...)
end

function substitute(ex::FlowExpression, this::TimeVariable, by::TimeExpression)
    FlowExpression(ex.fun, substitute(ex.t, this, by), substitute(ex.x, this, by), ex.dt_order,
        [substitute(x, this, by) for x in ex.d_args]...)
end

function substitute(ex::NonAutonomousFunctionExpression, this::TimeVariable, by::TimeExpression)
    NonAutonomousFunctionExpression(ex.fun, substitute(ex.t, this, by), substitute(ex.x, this, by), ex.dt_order,
        [substitute(x, this, by) for x in ex.d_args]...)
end

# substitute SpaceVariable by SpaceExpression

substitute(ex::SpaceVariable, this::SpaceVariable, by::SpaceExpression) = (ex==this ? by : ex)

function substitute(ex::SpaceLinearCombination, this::SpaceVariable, by::SpaceExpression)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

function substitute(ex::AutonomousFunctionExpression, this::SpaceVariable, by::SpaceExpression)
    AutonomousFunctionExpression(ex.fun, substitute(ex.x, this, by), [substitute(x, this, by) for x in ex.d_args]...)
end

function substitute(ex::FlowExpression, this::SpaceVariable, by::SpaceExpression)
    FlowExpression(ex.fun, ex.t, substitute(ex.x, this, by), ex.dt_order,
        [substitute(x, this, by) for x in ex.d_args]...)
end

function substitute(ex::NonAutonomousFunctionExpression, this::SpaceVariable, by::SpaceExpression)
    NonAutonomousFunctionExpression(ex.fun, ex.t, substitute(ex.x, this, by), ex.dt_order,
        [substitute(x, this, by) for x in ex.d_args]...)
end

# substitute AutonomousFunction by AutonomousFunction

substitute(ex::SpaceVariable, this::AutonomousFunction, by::AutonomousFunction) = ex

function substitute(ex::SpaceLinearCombination, this::AutonomousFunction, by::AutonomousFunction)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

function substitute(ex::AutonomousFunctionExpression, this::AutonomousFunction, by::AutonomousFunction)
    AutonomousFunctionExpression( ex.fun==this ? by : ex.fun, 
        substitute(ex.x, this, by), [substitute(x, this, by) for x in ex.d_args]...)
end

function substitute(ex::FlowExpression, this::AutonomousFunction, by::AutonomousFunction)
    FlowExpression(  ex.fun==this ? by : ex.fun,
        ex.t, substitute(ex.x, this, by), ex.dt_order, [substitute(x, this, by) for x in ex.d_args]...)
end

function substitute(ex::NonAutonomousFunctionExpression, this::AutonomousFunction, by::AutonomousFunction)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, substitute(ex.x, this, by), ex.dt_order, [substitute(x, this, by) for x in ex.d_args]...)
end

# substitute NonAutonomousFunction by NonAutonomousFunction

substitute(ex::SpaceVariable, this::NonAutonomousFunction, by::NonAutonomousFunction) = ex

function substitute(ex::SpaceLinearCombination, this::NonAutonomousFunction, by::NonAutonomousFunction)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

function substitute(ex::AutonomousFunctionExpression, this::NonAutonomousFunction, by::NonAutonomousFunction)
    AutonomousFunctionExpression( ex.fun, 
        substitute(ex.x, this, by), [substitute(x, this, by) for x in ex.d_args]...)
end

function substitute(ex::FlowExpression, this::NonAutonomousFunction, by::NonAutonomousFunction)
    FlowExpression(  ex.fun,
        ex.t, substitute(ex.x, this, by), ex.dt_order, [substitute(x, this, by) for x in ex.d_args]...)
end

function substitute(ex::NonAutonomousFunctionExpression, this::NonAutonomousFunction, by::NonAutonomousFunction)
    NonAutonomousFunctionExpression(  ex.fun==this ? by : ex.fun,
        ex.t, substitute(ex.x, this, by), ex.dt_order, [substitute(x, this, by) for x in ex.d_args]...)
end


#substitute AutonomousFunction by SpaceExpression

substitute(ex::SpaceVariable, this::AutonomousFunction, by::SpaceExpression, u::SpaceVariable) = ex

function substitute(ex::SpaceLinearCombination, this::AutonomousFunction, by::SpaceExpression, u::SpaceVariable)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(substitute(x, this, by, u), c) for (x, c) in ex.terms])
end

function substitute(ex::AutonomousFunctionExpression, this::AutonomousFunction, by::SpaceExpression, u::SpaceVariable)
    xx = substitute(ex.x, this, by, u)
    d_args = [substitute(x, this, by, u)  for x in ex.d_args]
    if this==ex.fun
        r = substitute(by, u, x_var)
        for d_arg in d_args
            r = differential(r, x_var, d_arg)
        end
        return substitute(r, x_var, xx)
    else
        return AutonomousFunctionExpression(ex.fun, xx, d_args...)
    end
end

function substitute(ex::FlowExpression, this::AutonomousFunction, by::SpaceExpression, u::SpaceVariable)
    FlowExpression(ex.fun, ex.t, substitute(ex.x, this, by, u), ex.dt_order,
        [substitute(x, this, by, u) for x in ex.d_args]...)
end

function substitute(ex::NonAutonomousFunctionExpression, this::AutonomousFunction, by::SpaceExpression, 
                    u::SpaceVariable)
    NonAutonomousFunctionExpression(ex.fun, ex.t, substitute(ex.x, this, by, u), ex.dt_order,
        [substitute(x, this, by, u) for x in ex.d_args]...)
end



# TODO: to be completed ( NonAutonomousFunction by expression, Function by zero

### commutator ########################################

#commutator(F::AutonomousFunction, G::AutonomousFunction, u::SpaceExpression) = F(u,G(u)) - G(u,F(u))
commutator(F::VectorFieldExpression, G::VectorFieldExpression, u::SpaceExpression) = F(u,G(u)) - G(u,F(u))

# double commutator [F,[G,H]] 
#function commutator(F::AutonomousFunction, G::AutonomousFunction, H::AutonomousFunction, u::SpaceExpression)
function commutator(F::VectorFieldExpression, G::VectorFieldExpression, H::VectorFieldExpression, u::SpaceExpression)
	v = SpaceVariable("v") 
        #X = AutonomousFunction("X")
        X = VectorFieldVariable("X")
        substitute(substitute(commutator(F,X,v), X , commutator(G,H,v), v), v, u)
end

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

function differential(ex::NonAutonomousFunctionExpression, with_respect_to::SpaceVariable,
                      applied_to::SpaceExpression)    
    terms = Array{SpaceExpression}(length(ex.d_args)+1)
    dx = differential(ex.x, with_respect_to, applied_to)   
    terms[1] = NonAutonomousFunctionExpression(ex.fun, ex.t, ex.x, ex.dt_order, vcat(ex.d_args, dx)...)
    for i=1:length(ex.d_args)
          d_args = copy(ex.d_args)
          d_args[i] = differential(ex.d_args[i], with_respect_to, applied_to)
          terms[i+1] = NonAutonomousFunctionExpression(ex.fun, ex.t, ex.x, ex.dt_order, d_args...)
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

function t_derivative(ex::FlowExpression, with_respect_to::TimeVariable; flag::Bool=false)    
    f = coefficient(ex.t, with_respect_to)
    dx = t_derivative(ex.x, with_respect_to, flag=flag)    
    terms = Array{SpaceExpression}(length(ex.d_args)+2)
    if f != 0
        if !flag
            t = AutonomousFunctionExpression(ex.fun, FlowExpression(ex.fun, ex.t, x_var, 0))
        else
            t = FlowExpression(ex.fun, ex.t, x_var, 0, AutonomousFunctionExpression(ex.fun, x_var))
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

function t_derivative(ex::NonAutonomousFunctionExpression, with_respect_to::TimeVariable; flag::Bool=false)    
    f = coefficient(ex.t, with_respect_to)
    dx = t_derivative(ex.x, with_respect_to, flag=flag)    
    terms = Array{SpaceExpression}(length(ex.d_args)+2)
    if f != 0
        #terms[1] = f*NonAutonomousFunctionExpression(ex.fun, ex.t, ex.x, ex.dt_order+1, vcat(ex.d_args, dx)...)
        terms[1] = f*NonAutonomousFunctionExpression(ex.fun, ex.t, ex.x, ex.dt_order+1, ex.d_args...)
    else
        terms[1] = x_zero
    end
    terms[2] = NonAutonomousFunctionExpression(ex.fun, ex.t, ex.x, ex.dt_order, vcat(ex.d_args, dx)...)
    for i=1:length(ex.d_args)
          d_args = copy(ex.d_args)
          d_args[i] = t_derivative(ex.d_args[i], with_respect_to, flag=flag)
          terms[i+2] = NonAutonomousFunctionExpression(ex.fun, ex.t, ex.x, ex.dt_order, d_args...)
    end      
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(x, 1) for x in terms])
end



### expand #############################################

expand(x::SpaceVariable) = x

function expand(ex::SpaceLinearCombination)   
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand(x), c) for (x, c) in ex.terms])
end

function _expander(ex::FunctionExpression, k::Integer)
    if k > length(ex.d_args)
        return ex
    end
    if !isa(ex.d_args[k], SpaceLinearCombination)
        return _expander(ex, k+1)
    end
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(_expander(
          isa(ex,AutonomousFunctionExpression)?
             AutonomousFunctionExpression(ex.fun, ex.x, vcat(ex.d_args[1:k-1], x, ex.d_args[k+1:end])...)
          :
          (isa(ex,NonAutonomousFunctionExpression)?
             NonAutonomousFunctionExpression(ex.fun, ex.t, ex.x, ex.dt_order, vcat(ex.d_args[1:k-1], 
                                             x, ex.d_args[k+1:end])...)
          :
             FlowExpression(ex.fun, ex.t, ex.x, ex.dt_order, vcat(ex.d_args[1:k-1], x, ex.d_args[k+1:end])...))
        ,1), c) for (x,c) in ex.d_args[k].terms])
end

function expand(ex::AutonomousFunctionExpression)   
    ex1 = AutonomousFunctionExpression(ex.fun, expand(ex.x), [expand(x) for x in ex.d_args]...)
    return _expander(ex1, 1)
end

function expand(ex::FlowExpression)   
    ex1 = FlowExpression(ex.fun, ex.t, expand(ex.x), ex.dt_order, [expand(x) for x in ex.d_args]...)
    return _expander(ex1, 1)
end

function expand(ex::NonAutonomousFunctionExpression)   
    ex1 = NonAutonomousFunctionExpression(ex.fun, ex.t, expand(ex.x), ex.dt_order, [expand(x) for x in ex.d_args]...)
    return _expander(ex1, 1)
end



###  expand_vector_field_expressions ##################

expand_vector_field_expressions(x::SpaceVariable) = x

function expand_vector_field_expressions(ex::SpaceLinearCombination)   
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand_vector_field_expressions(x), c) for (x, c) in ex.terms])
end

function expand_vector_field_expressions(ex::AutonomousFunctionExpression)   
    if ex.fun==op_zero
        return x_zero
    elseif isa(ex.fun, VectorFieldVariable)
        return   ex.fun(expand_vector_field_expressions(ex.x), [expand_vector_field_expressions(x) for x in ex.d_args]...)
    elseif isa(ex.fun, VectorFieldLinearCombination)
        return SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand_vector_field_expressions(t(ex.x, ex.d_args...)), c) 
                                      for (t, c) in ex.fun.terms])
    elseif isa(ex.fun, VectorFieldCommutator)
        A = ex.fun.A
        B = ex.fun.B
        u = ex.x 
        ex1 = expand_vector_field_expressions(A(u,B(u)) - B(u,A(u)))
        for arg in ex.d_args
            ex1 = differential(ex1, u, expand_vector_field_expressions(arg))
        end
        return ex1
    else
        @assert false "expected VectorFieldVariable, VectorFieldLinearCombination, or VectorFieldCommutator"
    end
end

function expand_vector_field_expressions(ex::FlowExpression)   
    FlowExpression(ex.fun, ex.t, expand_vector_field_expressions(ex.x), ex.dt_order, [expand_vector_field_expressions(x) for x in ex.d_args]...)
end

function expand_vector_field_expressions(ex::NonAutonomousFunctionExpression)   
    NonAutonomousFunctionExpression(ex.fun, ex.t, expand_vector_field_expressions(ex.x), ex.dt_order, [expand_vector_field_expressions(x) for x in ex.d_args]...)
end



###  expand_vector_field_exprssions ###################

expand_vector_field_expressions(x::SpaceVariable) = x

function expand_vector_field_expressions(ex::SpaceLinearCombination)   
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand(x), c) for (x, c) in ex.terms])
end



### reduce_order########################################

reduce_order(x::SpaceVariable) = x

function reduce_order(ex::SpaceLinearCombination)   
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(reduce_order(x), c) for (x, c) in ex.terms])
end

function reduce_order(ex::AutonomousFunctionExpression)   
    AutonomousFunctionExpression(ex.fun, reduce_order(ex.x), [reduce_order(x) for x in ex.d_args]...)
end

function reduce_order(ex::NonAutonomousFunctionExpression)   
    NonAutonomousFunctionExpression(ex.fun, ex.t, reduce_order(ex.x), ex.dt_order, 
                                   [reduce_order(x) for x in ex.d_args]...)
end

function _reduce_order_init(m::Integer)
    global ro_m = m 
    #global ro_F = AutonomousFunction("F")
    global ro_F = VectorFieldVariable("F")
    global ro_t = TimeVariable("t")
    global ro_u = SpaceVariable("u")
    global ro_vars = SpaceVariable[ SpaceVariable(string("v",k)) for k=1:m ]
    global ro_expr = Array{SpaceExpression}(m)
    for k=1:m
        ro_expr[k] = ( t_derivative(FlowExpression(ro_F, ro_t, ro_u, 0, ro_vars[1:k-1]...), ro_t, flag=false)
                      -t_derivative(FlowExpression(ro_F, ro_t, ro_u, 0, ro_vars[1:k-1]...), ro_t, flag=true)
                      +FlowExpression(ro_F, ro_t, ro_u, 0, AutonomousFunctionExpression(ro_F, ro_u), ro_vars[1:k-1]...) )
    end
end

function reduce_order(ex::FlowExpression)
    Fu = AutonomousFunctionExpression(ex.fun, ex.x)
    j = -1
    m = length(ex.d_args)
    for i = 1:m
        if ex.d_args[i]==Fu
            j=i
            break
        end
    end
    if j>=1 && m<=ro_m
        other_args = vcat(ex.d_args[1:j-1], ex.d_args[j+1:end])
        r = substitute(substitute(substitute(ro_expr[m], ro_F, ex.fun), ro_t, ex.t), ro_u, reduce_order(ex.x))
        for i = 1:m-1
            r = substitute(r, ro_vars[i], reduce_order(other_args[i]))
        end
        return reduce_order(r)
    else
        return FlowExpression(ex.fun, ex.t, reduce_order(ex.x), ex.dt_order, [reduce_order(x) for x in ex.d_args]...)
    end
end

### FE2DEF #####################################

FE2DEF(x::SpaceVariable) = x

function FE2DEF(ex::SpaceLinearCombination)   
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(FE2DEF(x), c) for (x, c) in ex.terms])
end

function FE2DEF(ex::AutonomousFunctionExpression) # TODO: check this !!!
    if isa(ex.x, FlowExpression) && ex.fun == ex.x.fun &&  length(ex.x.d_args)==0 && length(ex.d_args)==0
        xx = FE2DEF(ex.x.x)
        return FlowExpression(ex.fun, ex.x.t, xx, 0, AutonomousFunctionExpression(ex.fun, xx)) 
    elseif (isa(ex.x, FlowExpression) && ex.fun == ex.x.fun &&  length(ex.x.d_args)==0 && length(ex.d_args)==1
            && isa(ex.d_args[1], FlowExpression) && ex.fun == ex.d_args[1].fun && length(ex.d_args[1].d_args)==1
            && ex.x.x==ex.d_args[1].x )
        xx = FE2DEF(ex.x.x)
        v =  FE2DEF(ex.d_args[1].d_args[1])
        return (FlowExpression(ex.fun, ex.x.t, xx, 0,  AutonomousFunctionExpression(ex.fun, xx, v)) 
               +FlowExpression(ex.fun, ex.x.t, xx, 0,  AutonomousFunctionExpression(ex.fun, xx), v) ) 
    else
       return AutonomousFunctionExpression(ex.fun, FE2DEF(ex.x), [FE2DEF(x) for x in ex.d_args]...)
    end   
end

function FE2DEF(ex::FlowExpression)   
    FlowExpression(ex.fun, ex.t, FE2DEF(ex.x), ex.dt_order, [FE2DEF(x) for x in ex.d_args]...)
end

function FE2DEF(ex::NonAutonomousFunctionExpression)   
    NonAutonomousFunctionExpression(ex.fun, ex.t, FE2DEF(ex.x), ex.dt_order, [FE2DEF(x) for x in ex.d_args]...)
end

### DEF2FE #####################################

DEF2FE(x::SpaceVariable) = x

function DEF2FE(ex::SpaceLinearCombination)   
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(DEF2FE(x), c) for (x, c) in ex.terms])
end

function DEF2FE(ex::AutonomousFunctionExpression)   
    AutonomousFunctionExpression(ex.fun, DEF2FE(ex.x), [DEF2FE(x) for x in ex.d_args]...)
end

function DEF2FE(ex::FlowExpression)  # TODO: check this !!!
    if (length(ex.d_args)==1 && isa(ex.d_args[1], AutonomousFunctionExpression) &&   ex.fun==ex.d_args[1].fun 
       && ex.x==ex.d_args[1].x && length(ex.d_args[1].d_args)==0 )
        r = ex.d_args[1]
        return AutonomousFunctionExpression(r.fun, FlowExpression(ex.fun, ex.t, DEF2FE(ex.x), ex.dt_order))
    elseif (length(ex.d_args)==1 && isa(ex.d_args[1], AutonomousFunctionExpression) &&   ex.fun==ex.d_args[1].fun 
       && ex.x==ex.d_args[1].x && length(ex.d_args[1].d_args)==1 )
        r = ex.d_args[1]
        v = DEF2FE(ex.d_args[1].d_args[1])
        x = DEF2FE(ex.x)
        return (AutonomousFunctionExpression(r.fun, FlowExpression(ex.fun, ex.t, x, ex.dt_order), 
                                            FlowExpression(ex.fun, ex.t, x, ex.dt_order, v)) 
              - FlowExpression(ex.fun, ex.t, x, ex.dt_order, AutonomousFunctionExpression(r.fun, x) ,v))
    else
        return FlowExpression(ex.fun, ex.t, DEF2FE(ex.x), ex.dt_order, [DEF2FE(x) for x in ex.d_args]...)
    end
end

function DEF2FE(ex::NonAutonomousFunctionExpression)   
    NonAutonomousFunctionExpression(ex.fun, DEF2FE(ex.x), [DEF2FE(x) for x in ex.d_args]...)
end

