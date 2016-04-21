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

