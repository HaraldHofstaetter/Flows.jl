### substitute ########################################################

# substitute TimeExpression by TimeExpression

substitute(ex::LieDerivative, this::TimeExpression, by::TimeExpression) = ex

function substitute(ex::LieExponential, this::TimeExpression, by::TimeExpression)
    exp(substitute(ex.t, this, by), ex.DF, ex.label)
end

function substitute(ex::LieProduct, this::TimeExpression, by::TimeExpression)
    LieProduct(LieExpression[substitute(x, this, by) for x in ex.factors])
end

function substitute(ex::LieCommutator, this::TimeExpression, by::TimeExpression)
    LieCommutator(substitute(ex.A, this, by), substitute(ex.B, this, by))
end

function substitute(ex::LieLinearCombination, this::TimeExpression, by::TimeExpression)
    LieLinearCombination(Tuple{LieExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

# substitute VectorFieldExpression by VectorFieldExpression

function substitute(ex::LieDerivative, this::VectorFieldExpression, by::VectorFieldExpression)
    D(substitute(ex.F, this, by))
end

function substitute(ex::LieExponential, this::VectorFieldExpression, by::VectorFieldExpression)
    exp(ex.t, substitute(ex.DF, this, by), ex.label)
end

function substitute(ex::LieProduct, this::VectorFieldExpression, by::VectorFieldExpression)
    LieProduct(LieExpression[substitute(x, this, by) for x in ex.factors])
end

function substitute(ex::LieCommutator, this::VectorFieldExpression, by::VectorFieldExpression)
    LieCommutator(substitute(ex.A, this, by), substitute(ex.B, this, by))
end

function substitute(ex::LieLinearCombination, this::VectorFieldExpression, by::VectorFieldExpression)
    LieLinearCombination(Tuple{LieExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

# substitute LieExpression by LieExpression

substitute(ex::LieDerivative, this::LieExpression, by::LieExpression) = (ex == this ? by : ex)

function substitute(ex::LieExponential, this::LieExpression, by::LieExpression)
    ex == this ? by : exp(ex.t, substitute(ex.DF, this, by), ex.label)
end

function substitute(ex::LieProduct, this::LieExpression, by::LieExpression)
    ex == this ? by : LieProduct(LieExpression[substitute(x, this, by) for x in ex.factors])
end

function substitute(ex::LieCommutator, this::LieExpression, by::LieExpression)
    ex == this ? by : LieCommutator(substitute(ex.A, this, by), substitute(ex.B, this, by))
end

function substitute(ex::LieLinearCombination, this::LieExpression, by::LieExpression)
    ex == this ? by : LieLinearCombination(Tuple{LieExpression, Real}[(substitute(x, this, by), c) for (x, c) in ex.terms])
end

# LieExpressionToSpaceExpressionApplication

substitute(comb::LieExpressionToSpaceExpressionApplication, this::TimeExpression, by::TimeExpression) =
        apply(substitute(comb.lie_ex, this, by), substitute(comb.ex, this, by), comb.u)

substitute(comb::LieExpressionToSpaceExpressionApplication, this::SpaceExpression, by::SpaceExpression) =
        apply(comb.lie_ex, substitute(comb.ex, this, by), substitute(comb.u, this, by))

substitute(comb::LieExpressionToSpaceExpressionApplication, this::VectorFieldExpression, by::VectorFieldExpression) =
        apply(substitute(comb.lie_ex, this, by), substitute(comb.ex, this, by), comb.u)

substitute(comb::LieExpressionToSpaceExpressionApplication, this::LieExpression, by::LieExpression) =
        apply(substitute(comb.lie_ex, this, by), comb.ex, comb.u)


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

expand_lie_expressions(ex::LieDerivative; expand_commutators::Bool=true) = ex
expand_lie_expressions(ex::LieExponential; expand_commutators::Bool=true) = ex

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


function expand_lie_expressions(ex::LieCommutator; expand_commutators::Bool=true) 
    if expand_commutators
        ex1 = LieCommutator(expand_lie_expressions(ex.A, expand_commutators=expand_commutators), 
                            expand_lie_expressions(ex.B, expand_commutators=expand_commutators)) 
        return _commutator_expander(ex1, 1)
    else
        return ex
    end
end

function expand_lie_expressions(ex::LieLinearCombination; expand_commutators::Bool=true)
   LieLinearCombination(Tuple{LieExpression, Real}[(expand_lie_expressions(x, expand_commutators=expand_commutators), c) for (x, c) in ex.terms])
end

_expander_mul(a::LieExpression, b::LieExpression; expand_commutators::Bool=true) = a*b

function _expander_mul(a::LieLinearCombination, b::LieLinearCombination; expand_commutators::Bool=true)
    LieLinearCombination( 
        reshape(Tuple{LieExpression, Real}[ (expand_lie_expressions(x, 
            expand_commutators=expand_commutators)*expand_lie_expressions(y, expand_commutators=expand_commutators), c*d)   # CHECK !
            for (x,c) in a.terms, (y,d) in b.terms], 
            length(a.terms)*length(b.terms)) )
end   

function _expander_mul(a::LieExpression, b::LieLinearCombination; expand_commutators::Bool=true)
    LieLinearCombination( Tuple{LieExpression, Real}[ (_expander_mul(a,expand_lie_expressions(x, expand_commutators=expand_commutators)), c)  
    for (x,c) in b.terms] )  # CHECK !
end   

function _expander_mul(a::LieLinearCombination, b::LieExpression; expand_commutators::Bool=true)
    LieLinearCombination( Tuple{LieExpression, Real}[ (_expander_mul(expand_lie_expressions(x, expand_commutators=expand_commutators),b), c)  
    for (x,c) in a.terms] ) #CHECK !
end   

function expand_lie_expressions(ex::LieProduct; expand_commutators::Bool=true)
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

expand_lie_expressions(comb::LieExpressionToSpaceExpressionApplication; expand_commutators::Bool=true) = 
    apply(expand_lie_expressions(comb.lie_ex, expand_commutators=expand_commutators), comb.ex, comb.u)

#short forms:
expand(ex::LieExpression; expand_commutators::Bool=true) = 
    expand_lie_expressions(ex, expand_commutators=expand_commutators) 

expand(comb::LieExpressionToSpaceExpressionApplication; expand_commutators::Bool=true) = 
    expand_lie_expressions(comb, expand_commutators=expand_commutators)
#Note that 'expand' is also a method for SpaceExpression, so that expand could
#be called from parents of comb::LieExpressionToSpaceExpressionApplicatio <: SpaceExpression
#Then the lie pressions of comb will be expanded, which is probably not intended.
#But the convenience for calling a comb directly has priority.


## expand_lie_expressions for SpaceExpressions

expand_lie_expressions(ex::SpaceVariable; expand_commutators::Bool=true) = ex

function expand_lie_expressions(ex::SpaceLinearCombination; expand_commutators::Bool=true)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand_lie_expressions(x, expand_commutators=expand_commutators), c) 
    for (x, c) in ex.terms])
end

function expand_lie_expressions(ex::AutonomousFunctionExpression; expand_commutators::Bool=true)
    AutonomousFunctionExpression( ex.fun, 
        expand_lie_expressions(ex.x, expand_commutators=expand_commutators), [expand_lie_expressions(x, expand_commutators=expand_commutators) 
        for x in ex.d_args]...)
end

function expand_lie_expressions(ex::FlowExpression; expand_commutators::Bool=true)
    FlowExpression(  ex.fun,
        ex.t, expand_lie_expressions(ex.x, expand_commutators=expand_commutators), ex.dt_order, 
        [expand_lie_expressions(x, expand_commutators=expand_commutators) for x in ex.d_args]...)
end

function expand_lie_expressions(ex::NonAutonomousFunctionExpression; expand_commutators::Bool=true)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, expand_lie_expressions(ex.x, expand_commutators=expand_commutators), ex.dt_order, 
        [expand_lie_expressions(x, expand_commutators=expand_commutators) for x in ex.d_args]...)
end


### expand_lie_derivatives ##################

function expand_lie_derivatives(D::LieDerivative; expand_linear_combinations::Bool=true, expand_commutators::Bool=true)
    if expand_commutators && isa(D.F, VectorFieldCommutator)
        A = expand_lie_derivatives(LieDerivative(D.F.A), 
           expand_linear_combinations=expand_linear_combinations,
           expand_commutators=expand_commutators)
        B = expand_lie_derivatives(LieDerivative(D.F.B), 
           expand_linear_combinations=expand_linear_combinations, 
           expand_commutators=expand_commutators)
        return LieCommutator(B, A) # note the reversed order
    elseif expand_linear_combinations && isa(D.F, VectorFieldLinearCombination)
        return LieLinearCombination(Tuple{LieExpression, Real}[(expand_lie_derivatives(LieDerivative(x), 
           expand_linear_combinations=expand_linear_combinations, 
           expand_commutators=expand_commutators), c) for (x, c) in D.F.terms])
    else
        return D
    end
end

expand_lie_derivatives(ex::LieExponential; expand_linear_combinations::Bool=true, expand_commutators::Bool=true) = ex

function expand_lie_derivatives(ex::LieCommutator) 
    A = expand_lie_derivatives(ex.A, 
        expand_linear_combinations=expand_linear_combinations, 
        expand_commutators=expand_commutators)
    B = expand_lie_derivatives(ex.B, 
        expand_linear_combinations=expand_linear_combinations, 
        expand_commutators=expand_commutators)
    LieCommutator(A, B)
end    

function expand_lie_derivatives(ex::LieLinearCombination; expand_linear_combinations::Bool=true, expand_commutators::Bool=true)
   LieLinearCombination(Tuple{LieExpression, Real}[(expand_lie_derivatives(x, 
            expand_linear_combinations=expand_linear_combinations, 
            expand_commutators=expand_commutators), c) for (x, c) in ex.terms])
end

function expand_lie_derivatives(ex::LieProduct; expand_linear_combinations::Bool=true, expand_commutators::Bool=true)
    LieProduct(LieExpression[expand_lie_derivatives(x, 
        expand_linear_combinations=expand_linear_combinations, 
        expand_commutators=expand_commutators) for x in ex.factors])
end

expand_lie_derivatives(comb::LieExpressionToSpaceExpressionApplication; 
     expand_linear_combinations::Bool=true, expand_commutators::Bool=true) = 
     apply(expand_lie_derivatives(comb.lie_ex, 
         expand_linear_combinations=expand_linear_combinations, 
         expand_commutators=expand_commutators), comb.ex, comb.u)

expand_lie_derivatives(ex::LieExpression, expand_linear_combinations::Bool=true, expand_commutators::Bool=true) =
     expand_lie_derivatives(ex, 
         expand_linear_combinations=expand_linear_combinations, 
         expand_commutators=expand_commutators)

## expand_lie_derivatives for SpaceExpressions

expand_lie_derivatives(ex::SpaceVariable; expand_linear_combinations::Bool=true, expand_commutators::Bool=true) = ex

function expand_lie_derivatives(ex::SpaceLinearCombination; expand_linear_combinations::Bool=true, expand_commutators::Bool=true)
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(expand_lie_derivatives(x, 
        expand_linear_combinations=expand_linear_combinations, 
        expand_commutators=expand_commutators), c) for (x, c) in ex.terms])
end

function expand_lie_derivatives(ex::AutonomousFunctionExpression; expand_linear_combinations::Bool=true, expand_commutators::Bool=true)
    AutonomousFunctionExpression( ex.fun, 
        expand_lie_derivatives(ex.x, 
            expand_linear_combinations=expand_linear_combinations, 
            expand_commutators=expand_commutators), [expand_lie_derivatives(x, 
            expand_linear_combinations=expand_linear_combinations, 
            expand_commutators=expand_commutators) for x in ex.d_args]...)
end

function expand_lie_derivatives(ex::FlowExpression; expand_linear_combinations::Bool=true, expand_commutators::Bool=true)
    FlowExpression(  ex.fun,
        ex.t, expand_lie_derivatives(ex.x, 
            expand_linear_combinations=expand_linear_combinations, 
            expand_commutators=expand_commutators), ex.dt_order, [expand_lie_derivatives(x, 
            expand_linear_combinations=expand_linear_combinations, 
            expand_commutators=expand_commutators) for x in ex.d_args]...)
end

function expand_lie_derivatives(ex::NonAutonomousFunctionExpression; expand_linear_combinations::Bool=true, expand_commutators::Bool=true)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, expand_lie_derivatives(ex.x, 
            expand_linear_combinations=expand_linear_combinations, 
            expand_commutators=expand_commutators), ex.dt_order, [expand_lie_derivatives(x, 
            expand_linear_combinations=expand_linear_combinations, 
            expand_commutators=expand_commutators) for x in ex.d_args]...)
end


### merge_lie_derivatives ##################

merge_lie_derivatives(ex::LieDerivative; merge_linear_combinations::Bool=true, merge_commutators::Bool=true) = ex

merge_lie_derivatives(ex::LieExponential; merge_linear_combinations::Bool=true, merge_commutators::Bool=true) = ex

function merge_lie_derivatives(ex::LieCommutator; merge_linear_combinations::Bool=true, merge_commutators::Bool=true) 
    A = merge_lie_derivatives(ex.A, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators)
    B = merge_lie_derivatives(ex.B, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators)
    if merge_commutators && isa(A, LieDerivative) && isa(B, LieDerivative)
        return LieDerivative(VectorFieldCommutator(B.F, A.F))  # note the reversed order
    else
        return LieCommutator(A, B)
    end
end    

function merge_lie_derivatives(ex::LieLinearCombination; 
   merge_linear_combinations::Bool=true, merge_commutators::Bool=true)
       ex1 = LieLinearCombination(Tuple{LieExpression, Real}[(merge_lie_derivatives(x, 
             merge_linear_combinations=merge_linear_combinations, 
             merge_commutators=merge_commutators), c) for (x, c) in ex.terms])
       if merge_linear_combinations
           v_ex = op_zero
           r = lie_zero
           for (x,c) in ex1.terms
               if isa(x, LieDerivative)
                   v_ex = v_ex + c*x.F 
               else
                   r = r + c*x
               end
           end
           ex1 = r + D(v_ex)
       end
       return ex1
end

function merge_lie_derivatives(ex::LieProduct; 
    merge_linear_combinations::Bool=true, merge_commutators::Bool=true)
        LieProduct(LieExpression[merge_lie_derivatives(x, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators) for x in ex.factors])
end

merge_lie_derivatives(comb::LieExpressionToSpaceExpressionApplication; 
    merge_linear_combinations::Bool=true, merge_commutators::Bool=true) = 
        apply(merge_lie_derivatives(comb.lie_ex, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators), comb.ex, comb.u)

merge_lie_derivatives(ex::LieExpression;
    merge_linear_combinations::Bool=true, merge_commutators::Bool=true) = 
        merge_lie_derivatives(ex, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators)

## merge_lie_derivatives for SpaceExpressions

merge_lie_derivatives(ex::SpaceVariable; 
     merge_linear_combinations::Bool=true, merge_commutators::Bool=true) = ex

function merge_lie_derivatives(ex::SpaceLinearCombination; 
    merge_linear_combinations::Bool=true, merge_commutators::Bool=true)
        SpaceLinearCombination(Tuple{SpaceExpression, Real}[(merge_lie_derivatives(x, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators), c) for (x, c) in ex.terms])
end

function merge_lie_derivatives(ex::AutonomousFunctionExpression; 
    merge_linear_combinations::Bool=true, merge_commutators::Bool=true)
    AutonomousFunctionExpression( ex.fun, 
        merge_lie_derivatives(ex.x, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators), [merge_lie_derivatives(x, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators) for x in ex.d_args]...)
end

function merge_lie_derivatives(ex::FlowExpression; 
    merge_linear_combinations::Bool=true, merge_commutators::Bool=true)
    FlowExpression(  ex.fun,
        ex.t, merge_lie_derivatives(ex.x, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators), ex.dt_order, [merge_lie_derivatives(x, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators) for x in ex.d_args]...)
end

function merge_lie_derivatives(ex::NonAutonomousFunctionExpression; 
    merge_linear_combinations::Bool=true, merge_commutators::Bool=true)
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, merge_lie_derivatives(ex.x, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators), ex.dt_order, [merge_lie_derivatives(x, 
            merge_linear_combinations=merge_linear_combinations, 
            merge_commutators=merge_commutators) for x in ex.d_args]...)
end



### t_derivative ##################

t_derivative(ex::LieDerivative, t::TimeVariable; to_the_right::Array{Label, 1}=Label[R]) = lie_zero

function t_derivative(ex::LieLinearCombination, t::TimeVariable; to_the_right::Array{Label, 1}=Label[R])
   LieLinearCombination(Tuple{LieExpression, Real}[(t_derivative(x, t, to_the_right=to_the_right), c) for (x, c) in ex.terms])
end

function t_derivative(ex::LieCommutator, t::TimeVariable; to_the_right::Array{Label, 1}=Label[R])
    dA = t_derivative(ex.A, t, to_the_right=to_the_right)
    dB = t_derivative(ex.B, t, to_the_right=to_the_right)
    #LieCommutator(dA, ex.B) + LieCommutator(ex.A, dB)
    add_factorized(LieCommutator(dA, ex.B), LieCommutator(ex.A, dB))
end

function t_derivative(ex::LieExponential, t::TimeVariable; to_the_right::Array{Label, 1}=Label[R])
    c = coefficient(ex.t, t)
    ex.label in to_the_right ? c*(ex.DF*ex) : c*(ex*ex.DF)
end

function t_derivative(ex::LieProduct, t::TimeVariable; to_the_right::Array{Label, 1}=Label[R])
    if length(ex.factors) == 0
        return lie_zero
    elseif length(ex.factors) == 1
        return  t_derivative(ex.factors[1], t, to_the_right=to_the_right)
    else
        ex1 = ex.factors[1]
        ex2 = LieProduct(ex.factors[2:end])
        #return t_derivative(ex1,t)*ex2 + ex1*t_derivative(ex2,t)
        return add_factorized(t_derivative(ex1, t, to_the_right=to_the_right )*ex2,  
                               ex1*t_derivative(ex2, t, to_the_right=to_the_right))
    end
#    ex1 = lie_zero
#    for i in 1:length(ex.factors)
#        x = ex.factors[i]
#        if isa(x, LieExponential)
#            c = coefficient(x.t, t)
#            if to_the_right
#                ex1 = ex1 + c*LieProduct(vcat(ex.factors[1:i-1], x.DF ,ex.factors[i:end]))
#            else
#                ex1 = ex1 + c*LieProduct(vcat(ex.factors[1:i], x.DF, ex.factors[i+1:end]))
#            end
#        else
#            dx = t_derivative(x, t, to_the_right=to_the_right)
#            if dx != lie_zero
#                ex1 = ex1 + LieProduct(vcat(ex.factors[1:i-1], dx ,ex.factors[i+1:end]))
#            end
#        end
#    end
#    ex1
end

function t_derivative(comb::LieExpressionToSpaceExpressionApplication, t::TimeVariable) 
   apply(t_derivative(comb.lie_ex, t), comb.ex, comb.u) + apply(comb.lie_ex, t_derivative(comb.ex, t), comb.u)
end   


### normalize_lie_products ##################

function normalize_lie_products(ex::LieDerivative; to_the_right::Array{Label, 1}=Label[R])
    if isa(ex.F, VectorFieldLinearCombination) && length(ex.F.terms)==1
        c = ex.F.terms[1][2]
        term = ex.F.terms[1][1]
        if c==0 then
            return lie_zero
        else
            return c*D(term)
        end    
    else
        return ex
    end    
end

function normalize_lie_products(ex::LieLinearCombination; to_the_right::Array{Label, 1}=Label[R])
   LieLinearCombination(Tuple{LieExpression, Real}[(normalize_lie_products(x, to_the_right=to_the_right), c) for (x, c) in ex.terms])
end

function normalize_lie_products(ex::LieCommutator; to_the_right::Array{Label, 1}=Label[Rm])
    A = normalize_lie_products(ex.A, to_the_right=to_the_right)
    B = normalize_lie_products(ex.B, to_the_right=to_the_right)
    LieCommutator(A, B)
end

function normalize_lie_products(ex::LieExponential; to_the_right::Array{Label, 1}=Label[R])
    if ex.t == t_zero
        return lie_id
    elseif isa(ex.DF.F, VectorFieldLinearCombination) && length(ex.DF.F.terms)==1
        c = ex.DF.F.terms[1][2]
        term = ex.DF.F.terms[1][1]
        if c==0 then
            return lie_id
        else
            return exp(c*ex.t, term, ex.label)
        end    
    else
        return ex
    end    
end

function normalize_lie_products(comb::LieExpressionToSpaceExpressionApplication; to_the_right::Array{Label, 1}=Label[R]) 
   apply(normalize_lie_products(comb.lie_ex, to_the_right=to_the_right), comb.ex, comb.u) 
end   

function normalize_lie_products(ex::LieProduct; to_the_right::Array{Label, 1}=Label[R])
    ex1 = lie_id
    for factor in ex.factors
       ex1 = ex1 *  normalize_lie_products(factor, to_the_right=to_the_right)
    end
    c = 1
    if isa(ex1, LieLinearCombination) && length(ex1.terms)==1 
        if !isa(ex1.terms[1][1], LieProduct)
            return ex1
        end    
        c = ex1.terms[1][2]
        ex1 = ex1.terms[1][1]
    end 
      
    imax = length(ex1.factors)
    i = 1    
    ex2 = lie_id
    while i<=imax
        factor = ex1.factors[i]
        if isa(factor, LieDerivative) || isa(factor, LieExponential)
            n = 0
            label_left = default
            label_right = default
            t_left = t_zero
            t_right = t_zero
            F = false
            if isa(factor, LieDerivative)
                 F = factor.F
                 n = 1
            else
                F = factor.DF.F
                if factor.label in to_the_right
                    t_right = factor.t
                    label_right = factor.label
                else
                    t_left = factor.t
                    label_left = factor.label
                end    
            end
            i = i + 1

            while i<=imax && ((isa(ex1.factors[i], LieDerivative) && ex1.factors[i].F==F) 
                          || (isa(ex1.factors[i], LieExponential) &&  ex1.factors[i].DF.F==F) )

                factor = ex1.factors[i]
                if isa(factor, LieDerivative)
                    n = n + 1
                else
                    if factor.label in to_the_right
                        t_right = t_right + factor.t
                        label_right = factor.label # the most right occurence determines the label
                    else
                        t_left = t_left + factor.t
                    end    
                end
                
                i = i + 1
            end                          

            if t_left != t_zero 
                ex2 = ex2*exp(t_left, F, label_left)
            end
            ex2 = ex2 * D(F)^n
            if t_right != t_zero 
                ex2 = ex2*exp(t_right, F, label_right)
            end
        else
            ex2 = ex2*term
            i = i+1
        end    
    end
    return c*ex2
end

## normalize_lie_products for SpaceExpressions

normalize_lie_products(ex::SpaceVariable; to_the_right::Array{Label, 1}=Label[R]) = ex

function normalize_lie_products(ex::SpaceLinearCombination; to_the_right::Array{Label, 1}=Label[R])
    SpaceLinearCombination(Tuple{SpaceExpression, Real}[(normalize_lie_products(x, to_the_right=to_the_right), c) 
    for (x, c) in ex.terms])
end

function normalize_lie_products(ex::AutonomousFunctionExpression; to_the_right::Array{Label, 1}=Label[R])
    AutonomousFunctionExpression( ex.fun, 
        normalize_lie_products(ex.x, to_the_right=to_the_right), [normalize_lie_products(x, to_the_right=to_the_right) 
        for x in ex.d_args]...)
end

function normalize_lie_products(ex::FlowExpression; to_the_right::Array{Label, 1}=Label[R])
    FlowExpression(  ex.fun,
        ex.t, normalize_lie_products(ex.x, to_the_right=to_the_right), ex.dt_order, 
        [normalize_lie_products(x, to_the_right=to_the_right) for x in ex.d_args]...)
end

function normalize_lie_products(ex::NonAutonomousFunctionExpression; to_the_right::Array{Label, 1}=Label[R])
    NonAutonomousFunctionExpression(  ex.fun,
        ex.t, normalize_lie_products(ex.x, to_the_right=to_the_right), ex.dt_order, 
        [normalize_lie_products(x, to_the_right=to_the_right) for x in ex.d_args]...)
end



