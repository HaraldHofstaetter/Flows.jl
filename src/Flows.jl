VERSION >= v"0.4.5" && __precompile__()

module Flows

import 
    Base: (*), +, -, ^, string, show, write, writemime, expand, exp

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
export VectorFieldCommutator, op_zero, normalize, evaluate_vector_field_expressions
export @vector_fields


export Operator, OperatorExpression, OperatorLinearCombination
export D, LieExpression, LieDerivative, LieExponential, LieProduct
export LieCommutator, LieExSpaceExVarCombination, combine
export evaluate, evaluate_lie_expressions, expand_lie_expressions
export expand_commutators, expand_lie_commutators

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
    #global lie_zero = _register(LieLinearCombination(Tuple{SpaceExpression, Real}[], 0)) 
    _reduce_order_init(5)
end

end 
