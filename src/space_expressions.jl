immutable AutonomousFunction
   name::AbstractString
end

_str(f::AutonomousFunction; flat::Bool=false, latex::Bool=false) = f.name
string(f::AutonomousFunction) = f.name
show(io::IO, f::AutonomousFunction) = print(io, string(f))

###################################################################################################
# Essentially the same stuff as in time_expressions.jl with Time substituted by Space...

abstract SpaceExpression

immutable SpaceVariable <: SpaceExpression
   name::AbstractString
end 

_str(x::SpaceVariable; flat::Bool=false, latex::Bool=false) = x.name
string(x::SpaceVariable) = _str(x)
show(io::IO, x::SpaceVariable) = print(io, _str(x))


immutable SpaceLinearCombination <: SpaceExpression
   terms :: Array{Tuple{SpaceExpression, Real},1}
   function SpaceLinearCombination(terms :: Array{Tuple{SpaceExpression, Real},1}, dummy::Int)
       new(terms)
   end    
end

function SpaceLinearCombination(terms :: Array{Tuple{SpaceExpression, Real},1})
    d = Dict{SpaceExpression,Real}()
    for (x,c) in terms
        if isa(x, SpaceLinearCombination)
            # Nested SpaceLinearCombinations are expanded into parent SpaceLinearCombination
            for (x1, c1) in x.terms
                get!(d, x1, 0) 
                d[x1] += c*c1
            end
        else
            get!(d, x, 0) 
            d[x] += c
        end
    end
    for (key, val) in d
        # Terms with coefficient 0 are deleted
        if val == 0
            delete!(d, key)
        end
    end    
    if length(d) == 1 
        # Each linear combination consisting of only one term with coefficient 1 
        # is replaced by this term.
        # Note that this term has been already registered.
        for (key, val) in d
            if val==1
                return key
            end    
        end
    end    
    return _register(SpaceLinearCombination(Tuple{SpaceExpression, Real}[(key,val) for (key, val) in d], 0))
end

SpaceLinearCombination(x...) = SpaceLinearCombination(Tuple{SpaceExpression, Real}[(x[i],x[i+1]) for i=1:2:length(x) ])

+(a::SpaceExpression, b::SpaceExpression) = SpaceLinearCombination(a,1, b, 1)
-(a::SpaceExpression, b::SpaceExpression) = SpaceLinearCombination(a,1, b,-1)
-(a::SpaceExpression) = (-1)*a
*(f::Real, ex::SpaceExpression) = SpaceLinearCombination(ex,f)
*(f::Real, ex::SpaceLinearCombination) = SpaceLinearCombination( Tuple{SpaceExpression, Real}[ (x, f*c) for (x, c) in ex.terms ])
*(ex::SpaceExpression, f::Real) = f*ex

function _str(ex::SpaceLinearCombination; flat::Bool=false, latex::Bool=false) 
    if length(ex.terms) == 0 
        return "0"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            isa(x,SpaceLinearCombination)?"(":"", #TODO: which other types don't need parantheses?
            flat?_str_flat_arg_name(x):_str(x, latex=latex), 
            isa(x,SpaceLinearCombination)?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   

string(ex::SpaceLinearCombination) = _str(ex)
show(io::IO, ex::SpaceLinearCombination) = print(io, _str(ex))

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
        for d_arg in d_args            
            if d_arg==x_zero
                # If the AutonomousFunctionExpression is a derivative (i.e. a multilinear map),
                # and one argument of this map is 0 then the AutonomousFunctionExpression itself shall be zero
                return x_zero
            end
        end     
        _register(new(fun, t_zero, x, 0, SpaceExpression[d_arg for d_arg in d_args]))
    end    
end

#overloaded 'call' allows syntax 'F(x)' instead of 'AutonomousFunctionExpression(F, x)'
function call(F::AutonomousFunction, x::SpaceExpression, d_args...)
    AutonomousFunctionExpression(F::AutonomousFunction, x::SpaceExpression, d_args...) 
end

function _str(ex::AutonomousFunctionExpression; flat::Bool=false, latex::Bool=false)
    s = _str(ex.fun, latex=latex)
    n2 = length(ex.d_args)
    if latex
       if n2==1
          s = string(s, "'")
       elseif n2==2
          s = string(s, "''")
       elseif n2==3
          s = string(s, "'''")
       elseif n2>=4
          s = string(s, "^{(", n2, ")}")
       end
    elseif n2>0 
        s = string(s, "{", n2, "}")
    end
    s = string(s, latex?"(":"[",
        flat?_str_flat_arg_name(ex.x):_str(ex.x, latex=latex),
        latex?")":"]")
    if latex && n2==1
        x = ex.d_args[1]
        s = string(s, "\\cdot ", 
            isa(x, SpaceLinearCombination) ? "(":"",
            _str(x, latex=true),
            isa(x, SpaceLinearCombination) ? ")":"")
    elseif n2>0
        s = string(s, "(", 
        join([
            flat?_str_flat_arg_name(x):_str(x, latex=latex)
            for x in ex.d_args],","), 
        ")")
    end
    s
end

string(ex::AutonomousFunctionExpression) = _str(ex)
show(io::IO, ex::AutonomousFunctionExpression) = print(io, _str(ex))

immutable FlowExpression <: FunctionExpression
    fun::AutonomousFunction
    t::TimeExpression 
    x::SpaceExpression
    dt_order::Int 
    d_args::Array{SpaceExpression,1}
    function FlowExpression(fun::AutonomousFunction, 
                            t::TimeExpression, 
                            x::SpaceExpression, 
                            dt_order::Int,
                            d_args...)
        for d_arg in d_args            
            if d_arg==x_zero
                # If the FlowExpression is a derivative (i.e. a multilinear map),
                # and one argument of this map is 0 then the FlowExpression itself shall be zero
                return x_zero
            end
        end     
        _register(new(fun, t, x, dt_order, SpaceExpression[d_arg for d_arg in d_args]))
    end    
end

function E(fun::AutonomousFunction, t::TimeExpression, x::SpaceExpression, d_args...)  
    FlowExpression(fun, t, x, 0, d_args...) 
end

function _str(ex::FlowExpression; flat::Bool=false, latex::Bool=false)
    n1 = ex.dt_order
    n2 = length(ex.d_args)
    if latex
        s = ""
        if n1>=1
            s = string(s, "\\partial_{1}")
        end    
        if n1>=2
            s = string(s, "^{", n1, "}")
        end    
        if n2>=1
            s = string(s, "\\partial_{2}")
        end    
        if n2>=2
            s = string(s, "^{", n2, "}")
        end    
        s = string(s, "\\mathcal{E}_{", _str(ex.fun, latex=true), "}")
    else
        s = string("E_", _str(ex.fun))
        if n1>0 || n2>0 
            s = string(s, "{", n1, ",", n2, "}")
        end
    end
    s = string(s, latex?"(":"[",
        flat?_str_flat_arg_name(ex.t):_str(ex.t, latex=latex),
        ",",
        flat?_str_flat_arg_name(ex.x):_str(ex.x, latex=latex),
        latex?")":"]")
    if latex && n2==1
        x = ex.d_args[1]
        s = string(s, "\\cdot ", 
            isa(x, SpaceLinearCombination) ? "(":"",
            _str(x, latex=true),
            isa(x, SpaceLinearCombination) ? ")":"")
    elseif n2>0
        s = string(s, "(", 
        join([
            flat?_str_flat_arg_name(x):_str(x, latex=latex)
            for x in ex.d_args],","), 
        ")")
    end
    s
end

string(ex::FlowExpression) = _str(ex)
show(io::IO, ex::FlowExpression) = print(io, _str(ex))

writemime(io::IO, ::MIME"application/x-latex", ex::SpaceExpression) = write(io, "\$", _str(ex, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  ex::SpaceExpression) = write(io, "\$", _str(ex, latex=true), "\$")




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

function _get_register_key(ex::FlowExpression)
    string('E', _str_from_objref(ex.fun), ":", _str_from_objref(ex.x), ":", ex.dt_order, "|",
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


