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
# Essentially the same stuff as in time_expressions.jl with Time substituted by Space...

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


