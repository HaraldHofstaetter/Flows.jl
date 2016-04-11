module Flows

import 
    Base: (*), +, -, string, show

export TimeExpression, TimeVariable, TimeLinearCombination
export SpaceExpression, SpaceVariable, SpaceLinearCombination
export AutonomousFunction
export FunctionExpression, AutonomousFunctionExpression, FlowExpression
export coefficient, substitute, simplify
export print_time_expression_register, print_space_expression_register
export t_zero, x_zero

export @t_vars, @x_vars, @funs

export _str, _expand, _collect
export _register, _str_flat_arg_name, _get_register_key, _register
export _time_expression_index, _time_expression_register
export _space_expression_index, _space_expression_register


_str_from_objref(x) = hex(Int(pointer_from_objref(x)))

###################################################################################################

abstract TimeExpression

immutable TimeVariable <: TimeExpression
   name::AbstractString
   latex::AbstractString
end 

function TimeVariable(name::AbstractString; 
                      latex::AbstractString=name)
    TimeVariable(name, latex)
end   

string(t::TimeVariable) = t.name
show(io::IO, t::TimeVariable) = print(io, string(t))

immutable TimeLinearCombination <: TimeExpression
   terms :: Array{Tuple{TimeExpression, Real},1}
end

global t_zero = TimeLinearCombination([]) # "empty" TimeExpression

TimeLinearCombination(x...) = simplify(TimeLinearCombination([ (x[i],x[i+1]) for i=1:2:length(x) ]))

+(a::TimeExpression, b::TimeExpression) = TimeLinearCombination(a,1, b, 1)
-(a::TimeExpression, b::TimeExpression) = TimeLinearCombination(a,1, b,-1)
-(a::TimeExpression) = (-1)*a
*(f::Real, ex::TimeVariable) = TimeLinearCombination(ex,f)
*(f::Real, ex::TimeLinearCombination) = TimeLinearCombination( [ (x, f*c) for (x, c) in ex.terms ] )
*(ex::TimeExpression, f::Real) = f*ex

function _str(ex::TimeLinearCombination; flat::Bool=true) 
    if length(ex.terms) == 0 
        return "0"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)!=TimeVariable?"(":"", 
            flat?_str_flat_arg_name(x):string(x), 
            typeof(x)!=TimeVariable?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   

string(ex::TimeLinearCombination) = _str(ex, flat=false)
show(io::IO, ex::TimeLinearCombination) = print(io, string(ex))

function _expand(ex::TimeLinearCombination)
    TimeLinearCombination(
        vcat([ typeof(x)==TimeVariable ? (x, c) : (c*_expand(x)).terms for (x, c) in ex.terms ]...))
end

function _collect(ex::TimeLinearCombination)
    d = Dict(ex.terms)
    for key in keys(d)
        d[key] = 0
    end    
    for (key, val) in ex.terms
        d[key] += val
    end
    for (key, val) in d
        if val == 0
            delete!(d, key)
        end
    end    
    if length(d) == 1 
        #Each linear combination consisting of only one term with coefficient 1 
        #is replaced by this term.
        #Note that this term must have been already registered.
        for (key, val) in d
            if val==1
                return key
            end    
        end
    end    
    return _register(TimeLinearCombination([(key,val) for (key, val) in d]))
end

simplify(v::TimeVariable) = v
simplify(ex::TimeLinearCombination) = _collect(_expand(ex))

coefficient(ex::TimeVariable, v::TimeVariable) = (ex==v ? 1 : 0)

function coefficient(ex::TimeLinearCombination, v::TimeVariable)
    c = 0
    for (ex1, c1) in ex.terms
       c += c1*coefficient(ex1, v)
    end
    c
end 

substitute(ex::TimeVariable, this::TimeVariable, by::TimeExpression) = (ex==this ? by : ex)

function substitute(ex::TimeLinearCombination, this::TimeVariable, by::TimeExpression)
    TimeLinearCombination([(substitute(x, this, by), c) for (x, c) in ex.terms])
end

global _time_expression_index = Dict{TimeExpression,Int}()
global _time_expression_register = Dict{ASCIIString,Tuple{TimeExpression,Int}}()

function _get_register_key(ex::TimeLinearCombination)
    string('L', join([join([_str_from_objref(x), ':', c]) 
        for (x,c) in sort(ex.terms, 
            lt = (a,b) -> pointer_from_objref(a[1])<pointer_from_objref(b[1])) ],'|'))
end    

function _register(ex::TimeExpression)
    key = _get_register_key(ex)
    (ex,i) = get!(_time_expression_register, key) do
        i = length(_time_expression_index)
        _time_expression_index[ex] = i 
        (ex, i)
    end
    ex
end

_str_flat_arg_name(ex::TimeExpression) = string("#", _time_expression_index[ex])
_str_flat_arg_name(ex::TimeVariable) = ex.name

function print_time_expression_register()
    for (ex, i) in sort(collect(values(_time_expression_register)), lt = (a,b)-> a[2]<b[2])
        println("#", i, "\t" , _str(ex, flat=true))
     end
end

###################################################################################################

immutable AutonomousFunction
   name::AbstractString
   latex::AbstractString
end

function AutonomousFunction(name::AbstractString; 
                            latex::AbstractString=name)
    AutonomousFunction(name, latex)
end   

string(f::AutonomousFunction) = f.name
show(io::IO, f::AutonomousFunction) = print(io, string(f))

###################################################################################################
# Essentially the same stuff again with Time substituted by Space...

abstract SpaceExpression

immutable SpaceVariable <: SpaceExpression
   name::AbstractString
   latex::AbstractString
end 

function SpaceVariable(name::AbstractString; 
                      latex::AbstractString=name)
    SpaceVariable(name, latex)
end   

string(x::SpaceVariable) = x.name
show(io::IO, x::SpaceVariable) = print(io, string(x))


immutable SpaceLinearCombination <: SpaceExpression
   terms :: Array{Tuple{SpaceExpression, Real},1}
end

global x_zero = SpaceLinearCombination([]) # "empty" SpaceExpression

SpaceLinearCombination(x...) = simplify(SpaceLinearCombination([ (x[i],x[i+1]) for i=1:2:length(x) ]))

+(a::SpaceExpression, b::SpaceExpression) = SpaceLinearCombination(a,1, b, 1)
-(a::SpaceExpression, b::SpaceExpression) = SpaceLinearCombination(a,1, b,-1)
-(a::SpaceExpression) = (-1)*a
*(f::Real, ex::SpaceVariable) = SpaceLinearCombination(ex,f)
*(f::Real, ex::SpaceLinearCombination) = SpaceLinearCombination( [ (x, f*c) for (x, c) in ex.terms ] )
*(ex::SpaceExpression, f::Real) = f*ex

function _str(ex::SpaceLinearCombination; flat::Bool=true) 
    if length(ex.terms) == 0 
        return "0"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)!=SpaceVariable?"(":"", #TODO: which other types don't need parantheses?
            flat?_str_flat_arg_name(x):string(x), 
            typeof(x)!=SpaceVariable?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   

string(ex::SpaceLinearCombination) = _str(ex, flat=false)
show(io::IO, ex::SpaceLinearCombination) = print(io, string(ex))

abstract FunctionExpression <: SpaceExpression

immutable AutonomousFunctionExpression <: FunctionExpression
    fun::AutonomousFunction
    t::TimeExpression # = t_zero 
    x::SpaceExpression
    dt_order::Int # = 0
    d_args::Array{SpaceExpression,1}
    function AutonomousFunctionExpression(fun::AutonomousFunction, 
                                          x::SpaceExpression, 
                                          d_args...) 
        _register(new(fun, t_zero, x, 0, SpaceExpression[x for x in d_args]))
    end    
end

#overloaded 'call' allows syntax 'F(x)' instead of 'AutonomousFunctionExpression(F, x)'
function call(F::AutonomousFunction, x::SpaceExpression, d_args...)
    AutonomousFunctionExpression(F::AutonomousFunction, x::SpaceExpression, d_args...) 
end

function _str(ex::AutonomousFunctionExpression; flat::Bool=true)
    s = string(ex.fun)
    n2 = length(ex.d_args)
    if n2>0 
        s = string(s, "{", n2, "}")
    end
    s = string(s, "[",
        flat?_str_flat_arg_name(ex.x):string(ex.x),
        "]")
    if n2>0
        s = string(s, "(", 
        join([
            flat?_str_flat_arg_name(x):string(x)
            for x in ex.d_args],","), 
        ")")
    end
    s
end

string(ex::AutonomousFunctionExpression) = _str(ex, flat=false)
show(io::IO, ex::AutonomousFunctionExpression) = print(io, string(ex))

immutable FlowExpression <: FunctionExpression
    fun::AutonomousFunction
    t::TimeExpression 
    x::SpaceExpression
    dt_order::Int 
    d_args::Array{SpaceExpression,1}
end


substitute(ex::SpaceVariable, this::SpaceVariable, by::SpaceExpression) = (ex==this ? by : ex)
substitute(ex::SpaceVariable, this::TimeVariable, by::TimeExpression) = ex 

function substitute(ex::SpaceLinearCombination, this::SpaceVariable, by::SpaceExpression)
    SpaceLinearCombination([(substitute(x, this, by), c) for (x, c) in ex.terms])
end

function substitute(ex::SpaceLinearCombination, this::TimeVariable, by::TimeExpression)
    SpaceLinearCombination([(substitute(x, this, by), c) for (x, c) in ex.terms])
end


global _space_expression_index = Dict{SpaceExpression,Int}()
global _space_expression_register = Dict{ASCIIString,Tuple{SpaceExpression,Int}}()

function _get_register_key(ex::SpaceLinearCombination)
    string('L', join([join([_str_from_objref(x), ':', c]) 
        for (x,c) in sort(ex.terms, 
            lt = (a,b) -> pointer_from_objref(a[1])<pointer_from_objref(b[1])) ],'|'))
end    

function _get_register_key(ex::AutonomousFunctionExpression)
    string('A', _str_from_objref(ex.fun), ":", _str_from_objref(ex.x), "|",
        join([_str_from_objref(x)
            for x in sort(ex.d_args, 
            lt = (a,b) -> pointer_from_objref(a)<pointer_from_objref(b)) ],':'))
end    


function _register(ex::SpaceExpression)
    key = _get_register_key(ex)
    (ex,i) = get!(_space_expression_register, key) do
        i = length(_space_expression_index)
        _space_expression_index[ex] = i 
        (ex, i)
    end
    ex
end

_str_flat_arg_name(ex::SpaceExpression) = string("#", _space_expression_index[ex])
_str_flat_arg_name(ex::SpaceVariable) = ex.name

function print_space_expression_register()
    for (ex, i) in sort(collect(values(_space_expression_register)), lt = (a,b)-> a[2]<b[2])
        println("#", i, "\t" , _str(ex, flat=true))
     end
end

###################################################################################################


include("constructors.jl")

end #module Flows
