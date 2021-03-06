VERSION >= v"0.4.5" && __precompile__()

module Flows

import 
    Base: (*), +, -, ^, string, show, write, writemime, expand, apply, exp

export TimeExpression, TimeVariable, TimeLinearCombination
export SpaceExpression, SpaceVariable, SpaceLinearCombination
export FunctionObject, AutonomousFunction, NonAutonomousFunction
export FunctionExpression, AutonomousFunctionExpression, FlowExpression
export E, coefficient, substitute, commutator
export t_zero, x_zero
export differential, t_derivative, expand, reduce_order
export FE2DEF, DEF2FE

export @t_vars, @x_vars, @funs, @nonautonomous_funs

export VectorFieldExpression, VectorFieldVariable, VectorFieldLinearCombination
export VectorFieldCommutator, op_zero, normalize, resolve_vector_field_expressions

export D, LieExpression, LieDerivative, LieExponential, LieProduct
export LieCommutator, LieExpressionToSpaceExpressionApplication, apply
export evaluate_lie_expressions, evaluate, expand_lie_expressions
export expand_lie_commutators, expand_lie_derivatives
export merge_lie_derivatives, lie_zero, lie_id
export add_factorized, Label, default, L, R, show_exponential_labels
export normalize_lie_products

_str_from_objref(x) = hex(Int(pointer_from_objref(x)))

include("time_expressions.jl")
include("vector_field_expressions.jl")
include("space_expressions.jl")
include("constructors.jl")
include("library.jl")
include("lie_expressions.jl")
include("lie_library.jl")

function __init__()
    global t_zero = _register(TimeLinearCombination(Tuple{TimeExpression, Real}[],0)) 
    global op_zero = _register(VectorFieldLinearCombination(Tuple{VectorFieldExpression, Real}[], 0)) 
    global x_zero = _register(SpaceLinearCombination(Tuple{SpaceExpression, Real}[], 0)) 
    global lie_zero = _register(LieLinearCombination(Tuple{LieExpression, Real}[], 0)) 
    global lie_id =  _register(LieProduct(LieExpression[], 0))
    _reduce_order_init(5)
end

end 
