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
    evaluate_lie_expressions(expand_lie_commutators(lie_ex), ex, u)
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


### expand_lie_commutators ##############################

expand_lie_commutators(ex::LieDerivative) = ex
expand_lie_commutators(ex::LieExponential) = ex

function expand_lie_commutators(ex::LieCommutator) 
    A = expand_lie_commutators(ex.A)
    B = expand_lie_commutators(ex.B)
    A*B - B*A
end    

function expand_lie_commutators(ex::LieLinearCombination)
   LieLinearCombination(Tuple{LieExpression, Real}[(expand_lie_commutators(x), c) for (x, c) in ex.terms])
end

function expand_lie_commutators(ex::LieProduct)
    LieProduct(LieExpression[expand_lie_commutators(x) for x in ex.factors])
end

expand_lie_commutators(comb::LieExpressionToSpaceExpressionApplication) = apply(expand_lie_commutators(comb.lie_ex), comb.ex, comb.u)

## expand_lie_commutators for SpaceExpressions

expand_lie_commutators(ex::SpaceVariable) = ex

function expand_lie_commutators(ex::SpaceLinearCombination)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand_lie_commutators(x), c) for (x, c) in ex.terms])
end

function expand_lie_commutators(ex::AutonomousFunctionExpression)
    AutonomousFunctionExpression( ex.fun, 
        expand_lie_commutators(ex.x), [expand_lie_commutators(x) for x in ex.d_args]...)
end

function expand_lie_commutators(ex::FlowExpression)
    FlowExpression(  ex.fun,
        ex.t, expand_lie_commutators(ex.x), ex.dt_order, [expand_lie_commutators(x) for x in ex.d_args]...)
end

function expand_lie_commutators(ex::NonAutonomousFunctionExpression)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, expand_lie_commutators(ex.x), ex.dt_order, [expand_lie_commutators(x) for x in ex.d_args]...)
end


### expand_lie_expressions ############################################

expand_lie_expressions(ex::LieDerivative) = ex
expand_lie_expressions(ex::LieExponential) = ex
expand_lie_expressions(ex::LieCommutator) = commutator(ex.A, ex.B)

function expand_lie_expressions(ex::LieLinearCombination)
   LieLinearCombination(Tuple{LieExpression, Real}[(expand_lie_expressions(x), c) for (x, c) in ex.terms])
end

_expander_mul(a::LieExpression, b::LieExpression) = a*b

function _expander_mul(a::LieLinearCombination, b::LieLinearCombination)
    LieLinearCombination( 
        reshape(Tuple{LieExpression, Real}[ (expand_lie_expressions(x)*expand_lie_expressions(y), c*d)   # CHECK !
                for (x,c) in a.terms, (y,d) in b.terms], 
                length(a.terms)*length(b.terms)) )
end   

function _expander_mul(a::LieExpression, b::LieLinearCombination)
    LieLinearCombination( Tuple{LieExpression, Real}[ (_expander_mul(a,expand_lie_expressions(x)), c)  for (x,c) in b.terms] )  # CHECK !
end   

function _expander_mul(a::LieLinearCombination, b::LieExpression)
    LieLinearCombination( Tuple{LieExpression, Real}[ (_expander_mul(expand_lie_expressions(x),b), c)  for (x,c) in a.terms] ) #CHECK !
end   

function expand_lie_expressions(ex::LieProduct)
    if length(ex.factors)==0
        return ex
    end    
    ex1 = expand_lie_expressions(ex.factors[1])
    for i=2:length(ex.factors)
        ex1 = _expander_mul(ex1, ex.factors[i])
    end    
    return ex1 
end

expand_lie_expressions(comb::LieExpressionToSpaceExpressionApplication) = apply(expand_lie_expressions(comb.lie_ex), comb.ex, comb.u)

#short forms:
expand(ex::LieExpression) = expand_lie_expressions(ex) 

expand(comb::LieExpressionToSpaceExpressionApplication) = expand_lie_expressions(comb)
#Note that 'expand' is also a method for SpaceExpression, so that expand could
#be called from parents of comb::LieExpressionToSpaceExpressionApplicatio <: SpaceExpression
#Then the lie pressions of comb will be expanded, which is probably not intended.
#But the convenience for calling a comb directly has priority.


## expand_lie_expressions for SpaceExpressions

expand_lie_expressions(ex::SpaceVariable) = ex

function expand_lie_expressions(ex::SpaceLinearCombination)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand_lie_expressions(x), c) for (x, c) in ex.terms])
end

function expand_lie_expressions(ex::AutonomousFunctionExpression)
    AutonomousFunctionExpression( ex.fun, 
        expand_lie_expressions(ex.x), [expand_lie_expressions(x) for x in ex.d_args]...)
end

function expand_lie_expressions(ex::FlowExpression)
    FlowExpression(  ex.fun,
        ex.t, expand_lie_expressions(ex.x), ex.dt_order, [expand_lie_expressions(x) for x in ex.d_args]...)
end

function expand_lie_expressions(ex::NonAutonomousFunctionExpression)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, expand_lie_expressions(ex.x), ex.dt_order, [expand_lie_expressions(x) for x in ex.d_args]...)
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

