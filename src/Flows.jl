VERSION >= v"0.4.5" && __precompile__()

module Flows

import 
    Base: (*), +, -, string, show, write, writemime, expand, exp

export TimeExpression, TimeVariable, TimeLinearCombination
export SpaceExpression, SpaceVariable, SpaceLinearCombination
export FunctionObject, AutonomousFunction, NonAutonomousFunction
export FunctionExpression, AutonomousFunctionExpression, FlowExpression
export E, coefficient, substitute, commutator
export print_time_expression_register, print_space_expression_register
export t_zero, x_zero
export differential, t_derivative, expand, reduce_order
export FE2DEF, DEF2FE

export @t_vars, @x_vars, @funs, @nonautonomous_funs

export Operator, OperatorExpression, OperatorLinearCombination
export D, LieExpression, LieDerivative, LieExponential, LieMonomial
export LieExSpaceExVarCombination, combine, transform

_str_from_objref(x) = hex(Int(pointer_from_objref(x)))

include("time_expressions.jl")
include("space_expressions.jl")
include("constructors.jl")
include("library.jl")
include("lie_derivatives.jl")

function __init__()
    global t_zero = _register(TimeLinearCombination(Tuple{TimeExpression, Real}[],0)) 
    global x_zero = _register(SpaceLinearCombination(Tuple{SpaceExpression, Real}[], 0)) 
    _reduce_order_init(5)
end

end 
