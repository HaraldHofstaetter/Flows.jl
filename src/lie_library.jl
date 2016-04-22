### evaluate_lie_expressions ##########################################

function evaluate_lie_expressions(lie_ex::LieDerivative, ex::SpaceExpression, u::SpaceVariable)
    G = VectorFieldVariable("")
    ex1 = G(u, lie_ex.F(u))
    substitute(ex1, G, ex, u)    
end

function evaluate_lie_expressions(lie_ex::LieExponential, ex::SpaceExpression, u::SpaceVariable)
    G = VectorFieldVariable("")
    ex1 = G(E(lie_ex.DF.F, lie_ex.t, u))
    substitute(ex1, G, ex, u)    
end

function evaluate_lie_expressions(lie_ex::LieProduct, ex::SpaceExpression, u::SpaceVariable)
    ex1 = ex
    for x in reverse(lie_ex.factors)
        ex1 = evaluate_lie_expressions(x, ex1, u)
    end
    ex1
end

function evaluate_lie_expressions(lie_ex::LieCommutator, ex::SpaceExpression, u::SpaceVariable)
    evaluate_lie_expressions(evaluate_lie_commutators(lie_ex), ex, u)
end

function evaluate_lie_expressions(lie_ex::LieLinearCombination, ex::SpaceExpression, u::SpaceVariable)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(evaluate_lie_expressions(x, ex, u), c) for (x, c) in lie_ex.terms])
end

evaluate_lie_expressions(comb::LieExpressionToSpaceExpressionApplication) = evaluate_lie_expressions(comb.lie_ex, comb.ex, comb.u)
evaluate_lie_expressions(lie_ex::LieExpression, u::SpaceVariable) = evaluate_lie_expressions(lie_ex, u, u)

function  evaluate_lie_expressions(lie_ex::LieExpression, Fu::AutonomousFunctionExpression)
   @assert isa(Fu.x, SpaceVariable) string(Fu.fun, "(SpaceVariable) expected")
   evaluate_lie_expressions(lie_ex, Fu, Fu.x)
end

#short forms:
evaluate(comb::LieExpressionToSpaceExpressionApplication) = evaluate_lie_expressions(comb) 
evaluate(lie_ex::LieExpression, ex::SpaceExpression, u::SpaceVariable) = evaluate_lie_expressions(lie_ex, ex, u) 
evaluate(lie_ex::LieExpression, u::SpaceVariable) = evaluate_lie_expressions(lie_ex, u) 
evaluate(lie_ex::LieExpression, Fu::AutonomousFunctionExpression) = evaluate_lie_expressions(lie_ex, Fu)

## evaluate_lie_expressions for SpaceExpressions

evaluate_lie_expressions(ex::SpaceVariable) = ex

function evaluate_lie_expressions(ex::SpaceLinearCombination)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(evaluate_lie_expressions(x), c) for (x, c) in ex.terms])
end

function evaluate_lie_expressions(ex::AutonomousFunctionExpression)
    AutonomousFunctionExpression( ex.fun, 
        evaluate_lie_expressions(ex.x), [evaluate_lie_expressions(x) for x in ex.d_args]...)
end

function evaluate_lie_expressions(ex::FlowExpression)
    FlowExpression(  ex.fun,
        ex.t, evaluate_lie_expressions(ex.x), ex.dt_order, [evaluate_lie_expressions(x) for x in ex.d_args]...)
end

function evaluate_lie_expressions(ex::NonAutonomousFunctionExpression)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, evaluate_lie_expressions(ex.x), ex.dt_order, [evaluate_lie_expressions(x) for x in ex.d_args]...)
end


### evaluate_lie_commutators ##############################

evaluate_lie_commutators(ex::LieDerivative) = ex
evaluate_lie_commutators(ex::LieExponential) = ex

function evaluate_lie_commutators(ex::LieCommutator) 
    A = evaluate_lie_commutators(ex.A)
    B = evaluate_lie_commutators(ex.B)
    A*B - B*A
end    

function evaluate_lie_commutators(ex::LieLinearCombination)
   LieLinearCombination(Tuple{LieExpression, Real}[(evaluate_lie_commutators(x), c) for (x, c) in ex.terms])
end

function evaluate_lie_commutators(ex::LieProduct)
    LieProduct(LieExpression[evaluate_lie_commutators(x) for x in ex.factors])
end

evaluate_lie_commutators(comb::LieExpressionToSpaceExpressionApplication) = apply(evaluate_lie_commutators(comb.lie_ex), comb.ex, comb.u)

## evaluate_lie_commutators for SpaceExpressions

evaluate_lie_commutators(ex::SpaceVariable) = ex

function evaluate_lie_commutators(ex::SpaceLinearCombination)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(evaluate_lie_commutators(x), c) for (x, c) in ex.terms])
end

function evaluate_lie_commutators(ex::AutonomousFunctionExpression)
    AutonomousFunctionExpression( ex.fun, 
        evaluate_lie_commutators(ex.x), [evaluate_lie_commutators(x) for x in ex.d_args]...)
end

function evaluate_lie_commutators(ex::FlowExpression)
    FlowExpression(  ex.fun,
        ex.t, evaluate_lie_commutators(ex.x), ex.dt_order, [evaluate_lie_commutators(x) for x in ex.d_args]...)
end

function evaluate_lie_commutators(ex::NonAutonomousFunctionExpression)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, evaluate_lie_commutators(ex.x), ex.dt_order, [evaluate_lie_commutators(x) for x in ex.d_args]...)
end


### expand_lie_expressions ############################################

expand_lie_expressions(ex::LieDerivative; expand_commutators=true) = ex
expand_lie_expressions(ex::LieExponential; expand_commutators=true) = ex

function _commutator_expander(ex::LieExpression, k::Integer) #same technique as for VectorFieldCommutor
    if k > 2 || !isa(ex, LieCommutator)
        return ex
    end
    arg = (k==1?ex.A:ex.B)
    if !isa(arg, LieLinearCombination)
        return _commutator_expander(ex, k+1)
    end
    LieLinearCombination(Tuple{LieExpression, Real}[
        (_commutator_expander(k==1?LieCommutator(x, ex.B):LieCommutator(ex.A, x), 1 ), c) for (x,c) in arg.terms])
end


function expand_lie_expressions(ex::LieCommutator; expand_commutators=true) 
    if expand_commutators
        ex1 = LieCommutator(expand_lie_expressions(ex.A, expand_commutators=expand_commutators), 
                            expand_lie_expressions(ex.B, expand_commutators=expand_commutators)) 
        return _commutator_expander(ex1, 1)
    else
        return ex
    end
end

function expand_lie_expressions(ex::LieLinearCombination; expand_commutators=true)
   LieLinearCombination(Tuple{LieExpression, Real}[(expand_lie_expressions(x, expand_commutators=expand_commutators), c) for (x, c) in ex.terms])
end

_expander_mul(a::LieExpression, b::LieExpression; expand_commutators=true) = a*b

function _expander_mul(a::LieLinearCombination, b::LieLinearCombination; expand_commutators=true)
    LieLinearCombination( 
        reshape(Tuple{LieExpression, Real}[ (expand_lie_expressions(x, 
            expand_commutators=expand_commutators)*expand_lie_expressions(y, expand_commutators=expand_commutators), c*d)   # CHECK !
            for (x,c) in a.terms, (y,d) in b.terms], 
            length(a.terms)*length(b.terms)) )
end   

function _expander_mul(a::LieExpression, b::LieLinearCombination; expand_commutators=true)
    LieLinearCombination( Tuple{LieExpression, Real}[ (_expander_mul(a,expand_lie_expressions(x, expand_commutators=expand_commutators)), c)  
    for (x,c) in b.terms] )  # CHECK !
end   

function _expander_mul(a::LieLinearCombination, b::LieExpression; expand_commutators=true)
    LieLinearCombination( Tuple{LieExpression, Real}[ (_expander_mul(expand_lie_expressions(x, expand_commutators=expand_commutators),b), c)  
    for (x,c) in a.terms] ) #CHECK !
end   

function expand_lie_expressions(ex::LieProduct; expand_commutators=true)
    if length(ex.factors)==0
        return ex
    end    
    ex1 = expand_lie_expressions(ex.factors[1], expand_commutators=expand_commutators)
    for i=2:length(ex.factors)
        ex1 = _expander_mul(ex1, expand_lie_expressions(ex.factors[i], expand_commutators=expand_commutators),
              expand_commutators=expand_commutators)
    end    
    return ex1 
end

expand_lie_expressions(comb::LieExpressionToSpaceExpressionApplication; expand_commutators=true) = 
    apply(expand_lie_expressions(comb.lie_ex, expand_commutators=expand_commutators), comb.ex, comb.u)

#short forms:
expand(ex::LieExpression; expand_commutators=true) = 
    expand_lie_expressions(ex, expand_commutators=expand_commutators) 

expand(comb::LieExpressionToSpaceExpressionApplication; expand_commutators=true) = 
    expand_lie_expressions(comb, expand_commutators=expand_commutators)
#Note that 'expand' is also a method for SpaceExpression, so that expand could
#be called from parents of comb::LieExpressionToSpaceExpressionApplicatio <: SpaceExpression
#Then the lie pressions of comb will be expanded, which is probably not intended.
#But the convenience for calling a comb directly has priority.


## expand_lie_expressions for SpaceExpressions

expand_lie_expressions(ex::SpaceVariable; expand_commutators=true) = ex

function expand_lie_expressions(ex::SpaceLinearCombination; expand_commutators=true)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand_lie_expressions(x, expand_commutators=expand_commutators), c) 
    for (x, c) in ex.terms])
end

function expand_lie_expressions(ex::AutonomousFunctionExpression; expand_commutators=true)
    AutonomousFunctionExpression( ex.fun, 
        expand_lie_expressions(ex.x, expand_commutators=expand_commutators), [expand_lie_expressions(x, expand_commutators=expand_commutators) 
        for x in ex.d_args]...)
end

function expand_lie_expressions(ex::FlowExpression; expand_commutators=true)
    FlowExpression(  ex.fun,
        ex.t, expand_lie_expressions(ex.x, expand_commutators=expand_commutators), ex.dt_order, 
        [expand_lie_expressions(x, expand_commutators=expand_commutators) for x in ex.d_args]...)
end

function expand_lie_expressions(ex::NonAutonomousFunctionExpression; expand_commutators=true)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, expand_lie_expressions(ex.x, expand_commutators=expand_commutators), ex.dt_order, 
        [expand_lie_expressions(x, expand_commutators=expand_commutators) for x in ex.d_args]...)
end


### distribute_lie_derivatives ##################

function distribute_lie_derivatives(D::LieDerivative)
    if isa(D.F, VectorFieldCommutator)
        A = distribute_lie_derivatives(LieDerivative(D.F.A))
        B = distribute_lie_derivatives(LieDerivative(D.F.B))
        return LieCommutator(B, A) # note the reversed order
    elseif isa(D.F, VectorFieldLinearCombination)
        return LieLinearCombination(Tuple{LieExpression, Real}[(distribute_lie_derivatives(LieDerivative(x)), c) for (x, c) in D.F.terms])
    else
        return D
    end
end

distribute_lie_derivatives(ex::LieExponential) = ex

function distribute_lie_derivatives(ex::LieCommutator) 
    A = distribute_lie_derivatives(ex.A)
    B = distribute_lie_derivatives(ex.B)
    LieCommutator(A, B)
end    

function distribute_lie_derivatives(ex::LieLinearCombination)
   LieLinearCombination(Tuple{LieExpression, Real}[(distribute_lie_derivatives(x), c) for (x, c) in ex.terms])
end

function distribute_lie_derivatives(ex::LieProduct)
    LieProduct(LieExpression[distribute_lie_derivatives(x) for x in ex.factors])
end

distribute_lie_derivatives(comb::LieExpressionToSpaceExpressionApplication) = apply(distribute_lie_derivatives(comb.lie_ex), comb.ex, comb.u)
distribute_lie_derivatives(ex::LieExpression) = distribute_lie_derivatives(ex)

## distribute_lie_derivatives for SpaceExpressions

distribute_lie_derivatives(ex::SpaceVariable) = ex

function distribute_lie_derivatives(ex::SpaceLinearCombination)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(distribute_lie_derivatives(x), c) for (x, c) in ex.terms])
end

function distribute_lie_derivatives(ex::AutonomousFunctionExpression)
    AutonomousFunctionExpression( ex.fun, 
        distribute_lie_derivatives(ex.x), [distribute_lie_derivatives(x) for x in ex.d_args]...)
end

function distribute_lie_derivatives(ex::FlowExpression)
    FlowExpression(  ex.fun,
        ex.t, distribute_lie_derivatives(ex.x), ex.dt_order, [distribute_lie_derivatives(x) for x in ex.d_args]...)
end

function distribute_lie_derivatives(ex::NonAutonomousFunctionExpression)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, distribute_lie_derivatives(ex.x), ex.dt_order, [distribute_lie_derivatives(x) for x in ex.d_args]...)
end


### gather_lie_commutators ##################

gather_lie_commutators(ex::LieDerivative) = ex

gather_lie_commutators(ex::LieExponential) = ex

function gather_lie_commutators(ex::LieCommutator) 
    A = gather_lie_commutators(ex.A)
    B = gather_lie_commutators(ex.B)
    if isa(A, LieDerivative) && isa(B, LieDerivative)
        return LieDerivative(VectorFieldCommutator(B.F, A.F))
    else
        return LieCommutator(A, B)
    end
end    

function gather_lie_commutators(ex::LieLinearCombination)
   LieLinearCombination(Tuple{LieExpression, Real}[(gather_lie_commutators(x), c) for (x, c) in ex.terms])
end

function gather_lie_commutators(ex::LieProduct)
    LieProduct(LieExpression[gather_lie_commutators(x) for x in ex.factors])
end

gather_lie_commutators(comb::LieExpressionToSpaceExpressionApplication) = apply(gather_lie_commutators(comb.lie_ex), comb.ex, comb.u)
gather_lie_commutators(ex::LieExpression) = gather_lie_commutators(ex)

## gather_lie_commutators for SpaceExpressions

gather_lie_commutators(ex::SpaceVariable) = ex

function gather_lie_commutators(ex::SpaceLinearCombination)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(gather_lie_commutators(x), c) for (x, c) in ex.terms])
end

function gather_lie_commutators(ex::AutonomousFunctionExpression)
    AutonomousFunctionExpression( ex.fun, 
        gather_lie_commutators(ex.x), [gather_lie_commutators(x) for x in ex.d_args]...)
end

function gather_lie_commutators(ex::FlowExpression)
    FlowExpression(  ex.fun,
        ex.t, gather_lie_commutators(ex.x), ex.dt_order, [gather_lie_commutators(x) for x in ex.d_args]...)
end

function gather_lie_commutators(ex::NonAutonomousFunctionExpression)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, gather_lie_commutators(ex.x), ex.dt_order, [gather_lie_commutators(x) for x in ex.d_args]...)
end



### t_derivative ##################

t_derivative(ex::LieDerivative, t::TimeVariable; to_the_left::Bool=false) = lie_zero

function t_derivative(ex::LieLinearCombination, t::TimeVariable; to_the_left::Bool=false)
   LieLinearCombination(Tuple{LieExpression, Real}[(t_derivative(x, t, to_the_left=to_the_left), c) for (x, c) in ex.terms])
end

function t_derivative(ex::LieCommutator, t::TimeVariable; to_the_left::Bool=false)
    dA = t_derivative(ex.A, t, to_the_left=to_the_left)
    dB = t_derivative(ex.B, t, to_the_left=to_the_left)
    LieCommutator(dA, ex.B) + LieCommutator(ex.A, dB)
end

function t_derivative(ex::LieExponential, t::TimeVariable; to_the_left::Bool=false)
    c = coefficient(ex.t, t)
    to_the_left? c*(ex.DF*ex) : c*(ex*ex.DF)
end

function t_derivative(ex::LieProduct, t::TimeVariable; to_the_left::Bool=false)
    ex1 = lie_zero
    for i in 1:length(ex.factors)
        x = ex.factors[i]
        if isa(x, LieExponential)
            c = coefficient(x.t, t)
            if to_the_left
                ex1 = ex1 + c*LieProduct(vcat(ex.factors[1:i-1], x.DF ,ex.factors[i:end]))
            else
                ex1 = ex1 + c*LieProduct(vcat(ex.factors[1:i], x.DF, ex.factors[i+1:end]))
            end
        else
            dx = t_derivative(x, t, to_the_left=to_the_left)
            if dx != lie_zero
                ex1 = ex1 + LieProduct(vcat(ex.factors[1:i-1], dx ,ex.factors[i+1:end]))
            end
        end
    end
    ex1
end

function t_derivative(comb::LieExpressionToSpaceExpressionApplication, t::TimeVariable) 
   apply(t_derivative(comb.lie_ex, t), comb.ex, comb.u) + apply(comb.lie_ex, t_derivative(comb.ex, t), comb.u)
end   

