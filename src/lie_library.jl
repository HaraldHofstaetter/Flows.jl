### evaluate ##########################################

function evaluate(lie_ex::LieDerivative, ex::SpaceExpression, u::SpaceVariable)
    G = VectorFieldVariable("")
    if isa(lie_ex.F, VectorFieldLinearCombination)
        ex1 = G(u, SpaceLinearCombination( Tuple{SpaceExpression, Real}[ (F(u), c) 
                                           for (F, c) in lie_ex.F.terms ]) )
    else
        @assert isa(lie_ex.F, VectorFieldVariable) "internal error: VectorFieldVariable expected"
        ex1 = G(u, lie_ex.F(u))
    end    
    substitute(ex1, G, ex, u)    
end

function evaluate(lie_ex::LieExponential, ex::SpaceExpression, u::SpaceVariable)
    G = VectorFieldVariable("")
    if isa(lie_ex.DF.F, VectorFieldLinearCombination)
        if length(lie_ex.DF.F.terms)==1
            ex1 = G(E(lie_ex.DF.F.terms[1][1], lie_ex.t * lie_ex.DF.F.terms[1][2], u))
        else
            @assert false "Flow of VectorFieldLinearCombination not implemented"
        end
    else
        @assert isa(lie_ex.DF.F, VectorFieldVariable) "internal error: VectorFieldVariable expected"
        ex1 = G(E(lie_ex.DF.F, lie_ex.t, u))
    end    
    substitute(ex1, G, ex, u)    
end

function evaluate(lie_ex::LieProduct, ex::SpaceExpression, u::SpaceVariable)
    ex1 = ex
    for x in reverse(lie_ex.factors)
        ex1 = evaluate(x, ex1, u)
    end
    ex1
end

function evaluate(lie_ex::LieCommutator, ex::SpaceExpression, u::SpaceVariable)
    evaluate(expand_commutators(lie_ex), ex, u)
end

function evaluate(lie_ex::LieLinearCombination, ex::SpaceExpression, u::SpaceVariable)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(evaluate(x, ex, u), c) for (x, c) in lie_ex.terms])
end

evaluate(comb::LieExSpaceExVarCombination) = evaluate(comb.lie_ex, comb.ex, comb.u)

evaluate(lie_ex::LieExpression, u::SpaceVariable) = evaluate(lie_ex, u, u)

function  evaluate(lie_ex::LieExpression, Fu::AutonomousFunctionExpression)
   @assert isa(Fu.x, SpaceVariable) string(Fu.fun, "(SpaceVariable) expected")
   evaluate(lie_ex, Fu, Fu.x)
end

## evaluate_lie_expressions for SpaceExpressions

evaluate_lie_expressions(comb::LieExSpaceExVarCombination) = evaluate(comb)

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


### expand_commutators ##############################

expand_commutators(ex::LieDerivative) = ex
expand_commutators(ex::LieExponential) = ex

function expand_commutators(ex::LieCommutator) 
    A = expand_commutators(ex.A)
    B = expand_commutators(ex.B)
    A*B - B*A
end    

function expand_commutators(ex::LieLinearCombination)
   LieLinearCombination(Tuple{LieExpression, Real}[(expand_commutators(x), c) for (x, c) in ex.terms])
end

function expand_commutators(ex::LieProduct)
    LieProduct(LieExpression[expand_commutators(x) for x in ex.factors])
end

expand_commutators(comb::LieExSpaceExVarCombination) = combine(expand_commutators(comb.lie_ex), comb.ex, comb.u)
expand_lie_commutators(ex::LieExpression) = expand_commutators(ex)

## expand_lie_commutators for SpaceExpressions

expand_lie_commutators(comb::LieExSpaceExVarCombination) = expand_commutators(comb)
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


### expand ############################################

expand(ex::LieDerivative) = ex
expand(ex::LieExponential) = ex
expand(ex::LieCommutator) = commutator(ex.A, ex.B)

function expand(ex::LieLinearCombination)
   LieLinearCombination(Tuple{LieExpression, Real}[(expand(x), c) for (x, c) in ex.terms])
end

_expander_mul(a::LieExpression, b::LieExpression) = a*b

function _expander_mul(a::LieLinearCombination, b::LieLinearCombination)
    LieLinearCombination( 
        reshape(Tuple{LieExpression, Real}[ (x*y, c*d)   
                for (x,c) in a.terms, (y,d) in b.terms], 
                length(a.terms)*length(b.terms)) )
end   

function _expander_mul(a::LieExpression, b::LieLinearCombination)
    LieLinearCombination( Tuple{LieExpression, Real}[ (_expander_mul(a,x), c)  for (x,c) in b.terms] ) 
end   

function _expander_mul(a::LieLinearCombination, b::LieExpression)
    LieLinearCombination( Tuple{LieExpression, Real}[ (_expander_mul(x,b), c)  for (x,c) in a.terms] ) 
end   

function expand(ex::LieProduct)
    if length(ex.factors)==0
        return ex
    end    
    ex1 = expand(ex.factors[1])
    for i=2:length(ex.factors)
        ex1 = _expander_mul(ex1, ex.factors[i])
    end    
    return ex1 
end

expand(comb::LieExSpaceExVarCombination) = combine(expand(comb.lie_ex), comb.ex, comb.u)
expand_lie_expressions(ex::LieExpression) = expand(ex)

## expand_lie_expressions for SpaceExpressions

expand_lie_expressions(comb::LieExSpaceExVarCombination) = expand(comb)
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

