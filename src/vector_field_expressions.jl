abstract VectorFieldExpression # TODO: <: FunctionObject

immutable VectorFieldVariable <: VectorFieldExpression
    name::AbstractString
end

_str(t::VectorFieldVariable; flat::Bool=false, latex::Bool=false) = t.name


immutable VectorFieldLinearCombination <: VectorFieldExpression
    terms :: Array{Tuple{VectorFieldExpression, Real},1}
    function VectorFieldLinearCombination(terms::Array{Tuple{VectorFieldExpression, Real},1}, dummy::Int)
       # dummy only to make it distinguishable from the constructor below
       new(terms)
   end    
end

function VectorFieldLinearCombination(terms :: Array{Tuple{VectorFieldExpression, Real},1})
    d = Dict{VectorFieldExpression,Real}()
    for (x,c) in terms
        #@assert isa(c, Real) "Real expected"
        if isa(x, VectorFieldLinearCombination)
            # Nested VectorFieldLinearCombination are expanded into parent VectorFieldLinearCombination
            for (x1, c1) in x.terms
                #if isa(x1, VectorFieldCommutator)
                #    key = getkey(d, x1, false)
                #    if key==false
                #        key = getkey(d, VectorFieldCommutator(x1.B, x1.A), false)
                #        if key==false
                #            d[x1] = c*c1
                #        else
                #            d[key] -= c*c1
                #        end
                #    else
                #        d[x1] += c*c1
                #    end
                #else
                    get!(d, x1, 0)
                    d[x1] += c*c1
                #end
            end
        #elseif isa(x, VectorFieldCommutator)
        #    key = getkey(d, x, false)
        #    if key==false
        #        key = getkey(d, VectorFieldCommutator(x.B, x.A), false)
        #        if key==false
        #            d[x] = c
        #        else
        #            d[key] -= c
        #        end
        #    else
        #        d[x] += c
        #    end
        else
            #@assert isa(x, VectorFieldExpression) "VectorFieldExpression expected"
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
    return _register(VectorFieldLinearCombination(Tuple{VectorFieldExpression, Real}[(key,val) for (key, val) in d], 0))
end

VectorFieldLinearCombination(x...) = VectorFieldLinearCombination(Tuple{VectorFieldExpression, Real}[ (x[i],x[i+1]) for i=1:2:length(x) ])

+(a::VectorFieldExpression, b::VectorFieldExpression) = VectorFieldLinearCombination(a,1, b, 1)
-(a::VectorFieldExpression, b::VectorFieldExpression) = VectorFieldLinearCombination(a,1, b,-1)
-(a::VectorFieldExpression) = (-1)*a
+(a::VectorFieldExpression) = a
*(f::Real, ex::VectorFieldExpression) = VectorFieldLinearCombination(ex,f)
*(f::Real, ex::VectorFieldLinearCombination) = VectorFieldLinearCombination( Tuple{VectorFieldExpression, Real}[ (x, f*c) for (x, c) in ex.terms ] )
*(ex::VectorFieldExpression, f::Real) = f*ex


function _str(ex::VectorFieldLinearCombination; flat::Bool=false, latex::Bool=false) 
    if length(ex.terms) == 0 
        return latex?"0":"op_zero"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)==VectorFieldLinearCombination?"(":"", 
            flat?_str_flat_arg_name(x):_str(x, latex=latex), 
            typeof(x)==VectorFieldLinearCombination?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   


immutable VectorFieldCommutator <: VectorFieldExpression
    A::VectorFieldExpression
    B::VectorFieldExpression
    function VectorFieldCommutator (A::VectorFieldExpression, B::VectorFieldExpression)
        if A==B || A==op_zero || B==op_zero
            op_zero
        else
            _register(new(A,B))
        end    
    end        
end

commutator(A::VectorFieldExpression, B::VectorFieldExpression) = VectorFieldCommutator(A, B)

_str(ex::VectorFieldCommutator; flat::Bool=false, latex::Bool=false) = 
   string(latex?"[":"commutator(", 
          flat?_str_flat_arg_name(ex.A):_str(ex.A, latex=latex), ",", 
          flat?_str_flat_arg_name(ex.B):_str(ex.B, latex=latex), 
          latex?"]":")")

expand(x::VectorFieldVariable) = x

function expand(ex::VectorFieldLinearCombination)   
    VectorFieldLinearCombination(Tuple{VectorFieldExpression, Real}[(expand(x), c) for (x, c) in ex.terms])
end

function _expander(ex::VectorFieldExpression, k::Integer)
    if k > 2 || !isa(ex, VectorFieldCommutator)
        return ex
    end
    arg = (k==1?ex.A:ex.B)
    if !isa(arg, VectorFieldLinearCombination)
        return _expander(ex, k+1)
    end
    VectorFieldLinearCombination(Tuple{VectorFieldExpression, Real}[
        (_expander(k==1?VectorFieldCommutator(x, ex.B):VectorFieldCommutator(ex.A, x), 1 ), c) for (x,c) in arg.terms])
end

function expand(ex::VectorFieldCommutator)   
    ex1 = VectorFieldCommutator(expand(ex.A), expand(ex.B)) 
    return _expander(ex1, 1)
end

               
function _normalize(ex::VectorFieldExpression)
    if !isa(ex, VectorFieldCommutator) 
        return (ex, 1)
    end
    (A, sA) = _normalize(ex.A)
    (B, sB) = _normalize(ex.B) 
    s = sA*sB
    pointer_from_objref(A)<pointer_from_objref(B) ? 
        (VectorFieldCommutator(A,B), s) : (VectorFieldCommutator(B,A), -s)
end

function normalize(ex::VectorFieldExpression)
    ex1 = expand(ex)
    if isa(ex1, VectorFieldLinearCombination)
       terms = Tuple{Tuple{VectorFieldExpression, Int}, Real}[(_normalize(x),c) for (x,c) in ex1.terms ]
       return  VectorFieldLinearCombination(Tuple{VectorFieldExpression, Real}[ (x[1], x[2]*c) for (x,c) in terms])
    else
       (ex1, s) = _normalize(ex1)
       return s<1 ? -ex1 : ex1
    end
end


string(ex::VectorFieldExpression) = _str(ex)
show(io::IO, ex::VectorFieldExpression) = print(io, _str(ex))
writemime(io::IO, ::MIME"application/x-latex", ex::VectorFieldExpression) = write(io, "\$", _str(ex, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  ex::VectorFieldExpression) = write(io, "\$", _str(ex, latex=true), "\$")


global _vector_field_expression_index = Dict{VectorFieldExpression,Int}()
global _vector_field_expression_register = Dict{ASCIIString,Tuple{VectorFieldExpression,Int}}()

function _get_register_key(ex::VectorFieldLinearCombination)
    string('L', join([join([_str_from_objref(x), ':', c]) 
        for (x,c) in sort(ex.terms, 
            lt = (a,b) -> pointer_from_objref(a[1])<pointer_from_objref(b[1])) ],'|'))
end    

function _get_register_key(ex::VectorFieldCommutator)
    string('C', _str_from_objref(ex.A), ':', _str_from_objref(ex.B)) 
end    

function _register(ex::VectorFieldExpression)
    key = _get_register_key(ex)
    (ex,i) = get!(_vector_field_expression_register, key) do
        i = length(_vector_field_expression_index)
        _vector_field_expression_index[ex] = i 
        (ex, i)
    end
    ex
end

_str_flat_arg_name(ex::VectorFieldExpression) = string("#", _vector_field_expression_index[ex])
_str_flat_arg_name(ex::VectorFieldVariable) = ex.name

function print_vector_field_expression_register()
    for (ex, i) in sort(collect(values(_vector_field_expression_register)), lt = (a,b)-> a[2]<b[2])
        println("#", i, "\t" , _str(ex, flat=true))
     end
end


