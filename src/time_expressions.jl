abstract TimeExpression

immutable TimeVariable <: TimeExpression
   name::AbstractString
end 

_str(t::TimeVariable; flat::Bool=false, latex::Bool=false) = t.name
string(t::TimeVariable) = _str(t) 
show(io::IO, t::TimeVariable) = print(io, _str(t))

immutable TimeLinearCombination <: TimeExpression
   terms :: Array{Tuple{TimeExpression, Real},1}
   function TimeLinearCombination(terms :: Array{Tuple{TimeExpression, Real},1}, dummy::Int)
       # dummy only to make it distinguishable from the constructor below
       new(terms)
   end    
end

function TimeLinearCombination(terms :: Array{Tuple{TimeExpression, Real},1})
    d = Dict{TimeExpression,Real}()
    for (x,c) in terms
        #@assert isa(c, Real) "Real expected"
        if isa(x, TimeLinearCombination)
            # Nested TimeLinearCombinations are expanded into parent TimeLinearCombination
            for (x1, c1) in x.terms
                get!(d, x1, 0) 
                d[x1] += c*c1
            end
        else
            #@assert isa(x, TimeExpression) "TimeExpression expected"
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
    return _register(TimeLinearCombination(Tuple{TimeExpression, Real}[(key,val) for (key, val) in d], 0))
end


TimeLinearCombination(x...) = TimeLinearCombination(Tuple{TimeExpression, Real}[ (x[i],x[i+1]) for i=1:2:length(x) ])

+(a::TimeExpression, b::TimeExpression) = TimeLinearCombination(a,1, b, 1)
-(a::TimeExpression, b::TimeExpression) = TimeLinearCombination(a,1, b,-1)
-(a::TimeExpression) = (-1)*a
+(a::TimeExpression) = a
*(f::Real, ex::TimeVariable) = TimeLinearCombination(ex,f)
*(f::Real, ex::TimeLinearCombination) = TimeLinearCombination( Tuple{TimeExpression, Real}[ (x, f*c) for (x, c) in ex.terms ] )
*(ex::TimeExpression, f::Real) = f*ex

function _str(ex::TimeLinearCombination; flat::Bool=false, latex::Bool=false) 
    if length(ex.terms) == 0 
        return latex?"0":"t_zero"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)!=TimeVariable?"(":"", 
            flat?_str_flat_arg_name(x):_str(x, latex=latex), 
            typeof(x)!=TimeVariable?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   

string(ex::TimeLinearCombination) = _str(ex)
show(io::IO, ex::TimeLinearCombination) = print(io, _str(ex))

writemime(io::IO, ::MIME"application/x-latex", ex::TimeExpression) = write(io, "\$", _str(ex, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  ex::TimeExpression) = write(io, "\$", _str(ex, latex=true), "\$")


coefficient(ex::TimeVariable, v::TimeVariable) = (ex==v ? 1 : 0)

function coefficient(ex::TimeLinearCombination, v::TimeVariable)
    c = 0
    for (ex1, c1) in ex.terms
       c += c1*coefficient(ex1, v)
    end
    c
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

