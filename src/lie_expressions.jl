abstract LieExpression

immutable LieDerivative <: LieExpression
    F::VectorFieldExpression
    function LieDerivative(F::VectorFieldExpression)
        _register(new(F))
    end
end

D(F::VectorFieldExpression) = LieDerivative(F)

function _str(DF::LieDerivative; flat::Bool=false, latex::Bool=false) 
    if latex
        string("D_{", _str(DF.F, latex=true), "}")
    else   
        string("D(", _str(DF.F, flat= flat, latex=false), ")")
    end
end    


immutable Label
    name::AbstractString
end

global default = Label("")
global L = Label("L")
global R = Label("R")
global _show_exponential_labels = true

function show_exponential_labels(flag::Bool)
     global _show_exponential_labels = flag
end
show_exponential_labels()=show_exponential_labels(true)

immutable LieExponential <: LieExpression
    t ::TimeExpression
    DF::LieDerivative
    label::Label
    function LieExponential(t::TimeExpression, DF::LieDerivative, label::Label)
        _register(new(t, DF, label))
    end
end

exp(t::TimeExpression, DF::LieDerivative, label::Label) = LieExponential(t, DF, label)
exp(t::TimeExpression, DF::LieDerivative) = LieExponential(t, DF, default)
exp(t::TimeExpression, F::VectorFieldExpression, label::Label) =
    LieExponential(t, LieDerivative(F), label)
exp(t::TimeExpression, F::VectorFieldExpression) =
    LieExponential(t, LieDerivative(F), default)

function _str(E::LieExponential; flat::Bool=false, latex::Bool=false) 
    if latex
        string(_show_exponential_labels&&E.label!=default?"\\underbrace{":"",
               "\\mathrm{e}^{",  
               typeof(E.t)==TimeLinearCombination&&length(E.t.terms)>1?"(":"", 
              _str(E.t, latex=true),
               typeof(E.t)==TimeLinearCombination&&length(E.t.terms)>1?")":"", 
              _str(E.DF, latex=true), "}",
              _show_exponential_labels&&E.label!=default?string("}_\\mathrm{",E.label.name,"}"):"")
    else   
        string("exp(", _str(E.t, flat=flat, latex=false), ",",
               _str(E.DF, flat=flat, latex=false), 
               E.label==default?"":string(",",E.label.name),
               ")")
    end
end    

immutable LieCommutator <: LieExpression
    A::LieExpression
    B::LieExpression
    function LieCommutator(A::LieExpression, B::LieExpression)
        if A==B || A==lie_zero || B==lie_zero
            lie_zero
        else
            _register(new(A,B))
        end    
    end        
end

commutator(A::LieExpression, B::LieExpression) = LieCommutator(A, B)

_str(ex::LieCommutator; flat::Bool=false, latex::Bool=false) = 
   string(latex?"[":"commutator(", 
          flat?_str_flat_arg_name(ex.A):_str(ex.A, latex=latex), ",", 
          flat?_str_flat_arg_name(ex.B):_str(ex.B, latex=latex), 
          latex?"]":")")


immutable LieLinearCombination <: LieExpression 
    terms :: Array{Tuple{LieExpression, Real},1}
    function LieLinearCombination(terms
        :: Array{Tuple{LieExpression, Real},1}, dummy::Int)
       # dummy only to make it distinguishable from the constructor below
       new(terms)
   end    
end

function LieLinearCombination(terms :: Array{Tuple{LieExpression, Real},1})
    d = Dict{LieExpression,Real}()
    for (x,c) in terms
        #@assert isa(c, Real) "Real expected"
        if isa(x, LieLinearCombination)
            # Nested LieLinearCombination are expanded into parent LieLinearCombination
            for (x1, c1) in x.terms
                get!(d, x1, 0) 
                d[x1] += c*c1
            end
        else
            #@assert isa(x, LieExpression) "LieExpression expected"
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
    return _register(LieLinearCombination(Tuple{LieExpression, Real}[(key,val) for (key, val) in d], 0))
end

LieLinearCombination(x...) = LieLinearCombination(Tuple{LieExpression, Real}[ (x[i],x[i+1]) for i=1:2:length(x) ])

+(a::LieExpression, b::LieExpression) = LieLinearCombination(a,1, b, 1)
-(a::LieExpression, b::LieExpression) = LieLinearCombination(a,1, b,-1)
-(a::LieExpression) = (-1)*a
+(a::LieExpression) = a
*(f::Real, ex::LieExpression) = LieLinearCombination(ex,f)
*(f::Real, ex::LieLinearCombination) = LieLinearCombination( Tuple{LieExpression, Real}[ (x, f*c) for (x, c) in ex.terms ] )
*(ex::LieExpression, f::Real) = f*ex

##The following code would expand Products of LinearCombinations
##immediately, which is not intended and thus commented out.
#function *(a::LieLinearCombination, b::LieLinearCombination)
#    LieLinearCombination( 
#        reshape(Tuple{LieExpression, Real}[ (x*y, c*d)   
#                for (x,c) in a.terms, (y,d) in b.terms], 
#                length(a.terms)*length(b.terms)) )
#end   
#
#function *(a::Union{LieDerivative,LieExponential,LieProduct}, b::LieLinearCombination)
#    LieLinearCombination( Tuple{LieExpression, Real}[ (a*x, c)  for (x,c) in b.terms] ) 
#end   
#
#function *(a::LieLinearCombination, b::Union{LieDerivative,LieExponential,LieProduct})
#    LieLinearCombination( Tuple{LieExpression, Real}[ (x*b, c)  for (x,c) in a.terms] ) 
#end   


function _str(ex::LieLinearCombination; flat::Bool=false, latex::Bool=false) 
    if length(ex.terms) == 0 
        return latex?"0":"lie_zero"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)==LieLinearCombination?"(":"", 
            flat?_str_flat_arg_name(x):_str(x, latex=latex), 
            typeof(x)==LieLinearCombination?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   



immutable LieProduct <: LieExpression
    factors::Array{LieExpression,1}
    function LieProduct(factors::Array{LieExpression,1}, dummy::Int)
        new(factors)
    end    
end

function LieProduct(factors::Array{LieExpression,1})
    if length(factors)==0
        return lie_id
    elseif length(factors)==1
        return factors[1]
    else    
        return _register(LieProduct(factors, 0))
    end    
end


function _str(p; flat::Bool=false, latex::Bool=false) 
    if length(p.factors)==0
        return latex?"\\mathrm{Id}":"lie_id"
    end
    join([string(
        typeof(x)==LieLinearCombination?"(":"", 
        flat?_str_flat_arg_name(x):_str(x, latex=latex), 
        typeof(x)==LieLinearCombination?")":"")
    for x in p.factors], latex?"":"*")
end   

LieExprButNotProduct = Union{LieDerivative,LieExponential,LieLinearCombination, LieCommutator}
#With this trick the construction of Products of Products using only the operator * 
#(not the default constructor) is not possible (i.e. products are flattened).

#apply simple rules 0*a = a*0 = 0, 1*a = a*1 = a
*(a::LieExprButNotProduct, b::LieExprButNotProduct) = 
    a==lie_zero||b==lie_zero ? lie_zero : LieProduct(LieExpression[a, b])
*(a::LieProduct, b::LieProduct) = 
    a==lie_id ? b : ( b==lie_id ? a :LieProduct(vcat(a.factors, b.factors)) )
*(a::LieExprButNotProduct, b::LieProduct) =
    a==lie_zero ? lie_zero : ( b==lie_id ? a : LieProduct(vcat(a, b.factors)) )
*(a::LieProduct, b::LieExprButNotProduct) = 
    b==lie_zero ? lie_zero : ( a==lie_id ? b : LieProduct(vcat(a.factors, b)) )

function *(a::LieLinearCombination, b::LieLinearCombination)
   if length(a.terms)==1 && length(b.terms)==1
       return (a.terms[1][2]*b.terms[1][2])*(a.terms[1][1]*b.terms[1][1])
   else
       return a==lie_zero||b==lie_zero ? lie_zero : LieProduct(LieExpression[a, b])
   end
end

function *(a::LieLinearCombination, b::LieProduct)
    if length(a.terms)==1 
        return b==lie_id ? a : a.terms[1][2]*LieProduct(vcat(a.terms[1][1], b.factors)) 
    else
        return a==lie_zero ? lie_zero : LieProduct(vcat(a, b.factors)) 
    end
end

function *(a::LieLinearCombination, b::LieExprButNotProduct)
    if length(a.terms)==1 
        return a.terms[1][2]*(a.terms[1][1]*b)
    else
        return a==lie_zero||b==lie_zero ? lie_zero : LieProduct(LieExpression[a, b])
    end
end

function *(a::LieProduct, b::LieLinearCombination)
    if length(b.terms)==1 
        return a==lie_id ? b : b.terms[1][2]*LieProduct(vcat(a.factors, b.terms[1][1])) 
    else
        return b==lie_zero ? lie_zero : LieProduct(vcat(a.factors, b)) 
    end
end

function *(a::LieExprButNotProduct, b::LieLinearCombination)
    if length(b.terms)==1 
        return b.terms[1][2]*(a*b.terms[1][1])
    else
        return a==lie_zero||b==lie_zero ? lie_zero : LieProduct(LieExpression[a, b])
    end
end


^(a::LieDerivative, p::Integer)= LieProduct(LieExpression[a for i=1:p])    


add_factorized(a::LieExpression, b::LieExpression; ca::Real=1, cb::Real=1) = ca*a + cb*b

function add_factorized(a::LieProduct, b::LieProduct; ca::Real=1, cb::Real=1) 
    if a==b
       return 2*a
    end   
    h = min(length(a.factors), length(b.factors))
    i = 0
    while i<h && a.factors[i+1]==b.factors[i+1]
        i+=1
    end
    j = 0
    while j<h && a.factors[end-j]==b.factors[end-j]
        j+=1
    end    
    return LieProduct(vcat(a.factors[1:i], 
               ca*LieProduct(a.factors[i+1:end-j]) + cb*LieProduct(b.factors[i+1:end-j]),
               a.factors[end-j+1:end]))
end

function add_factorized(a::LieExpression, b::LieProduct; ca::Real=1, cb::Real=1) 
   if length(b.factors)>=1 && b.factors[1]==a
       return a*(ca*lie_id + cb*LieProduct(b.factors[2:end]))
   elseif length(b.factors)>=2 && b.factors[end]==a
       return (ca*lie_id+cb*LieProduct(b.factors[1:end-1]))*a
   else
       return ca*a+cb*b
   end
end   

add_factorized(a::LieProduct, b::LieExpression; ca::Real=1, cb::Real=1) = 
    add_factorized(b, a, ca=cb, cb=ca) 

function add_factorized(a::LieLinearCombination, b::LieLinearCombination; ca::Real=1, cb::Real=1) 
    if length(a.terms)==1 && length(b.terms)==1
        return add_factorized(a.terms[1][1], b.terms[1][1], ca=ca*a.terms[1][2], cb=cb*b.terms[1][2])
    else
        return ca*a+cb*b
    end    
end

#disambigouize 
function add_factorized(a::LieLinearCombination, b::LieProduct; ca::Real=1, cb::Real=1) 
    if length(a.terms)==1 
        return add_factorized(a.terms[1][1], b, ca=ca*a.terms[1][2], cb=cb)
    elseif length(b.factors)>=1 && b.factors[1]==a
        return a*(ca*lie_id + cb*LieProduct(b.factors[2:end]))
    elseif length(b.factors)>=2 && b.factors[end]==a
        return (ca*lie_id+cb*LieProduct(b.factors[1:end-1]))*a
    else
        return ca*a+cb*b
    end
end   

add_factorized(a::LieProduct, b::LieLinearCombination; ca::Real=1, cb::Real=1) =
    add_factorized(b, a, ca=cb, cb=ca) 

function add_factorized(a::LieLinearCombination, b::LieExpression; ca::Real=1, cb::Real=1) 
    if length(a.terms)==1 
        return add_factorized(a.terms[1][1], b, ca=ca*a.terms[1][2], cb=cb)
    else
        return ca*a+cb*b
    end    
end

add_factorized(a::LieExpression, b::LieLinearCombination; ca::Real=1, cb::Real=1) =
    add_factorized(b, a, ca=cb, cb=ca) 




string(ex::LieExpression) = _str(ex)
show(io::IO, ex::LieExpression) = print(io, _str(ex))

writemime(io::IO, ::MIME"application/x-latex", ex::LieExpression) = write(io, "\$", _str(ex, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  ex::LieExpression) = write(io, "\$", _str(ex, latex=true), "\$")




immutable LieExpressionToSpaceExpressionApplication <: SpaceExpression
    lie_ex::LieExpression
    ex::SpaceExpression
    u::SpaceVariable
end
 
apply(lie_ex::LieExpression, ex::SpaceExpression, u::SpaceVariable) =
   LieExpressionToSpaceExpressionApplication(lie_ex::LieExpression, ex::SpaceExpression, u::SpaceVariable)

apply(lie_ex::LieExpression, u::SpaceVariable) = apply(lie_ex, u, u)

function  apply(lie_ex::LieExpression, Fu::AutonomousFunctionExpression)
   @assert isa(Fu.x, SpaceVariable) string(Fu.fun, "(SpaceVariable) expected")
   apply(lie_ex, Fu, Fu.x)
end

*(lie_ex::LieExpression, u::SpaceVariable) = apply(lie_ex, u)
*(lie_ex::LieExpression, Fu::AutonomousFunctionExpression) = apply(lie_ex, Fu)


function _str(comb::LieExpressionToSpaceExpressionApplication; flat::Bool=false, latex::Bool=false, arg::SpaceVariable=_no_x_var) 
    if latex
        s = string(
            isa(comb.lie_ex, LieLinearCombination)?"(":"",
                _str(comb.lie_ex, flat=flat, latex=true),
            isa(comb.lie_ex, LieLinearCombination)?")":"",
            )
            if comb.ex==comb.u
                s = string(s,"\\mathrm{Id}(",_str(comb.u, flat=flat, latex=true), ")")
            elseif isa(comb.ex, AutonomousFunctionExpression) && comb.ex.x==comb.u
                s = string(s, _str(comb.ex, flat=flat, latex=true))
            else
                s = string(s,"[", _str(comb.ex, flat=flat, latex=true, arg=comb.u), 
                           "](", _str(comb.u, flat=flat, latex=true) , ")")
            end
    else
        s = string("apply(",
            _str(comb.lie_ex, flat=flat, latex=false),
            ",",
            _str(comb.ex, flat=flat, latex=false),
            ",",
            _str(comb.u, flat=flat, latex=false),
            ")")
    end
    s
end

string(comb::LieExpressionToSpaceExpressionApplication) = _str(comb)
show(io::IO, comb::LieExpressionToSpaceExpressionApplication) = print(io, _str(comb))
writemime(io::IO, ::MIME"application/x-latex", comb::LieExpressionToSpaceExpressionApplication) = write(io, "\$", _str(comb, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  comb::LieExpressionToSpaceExpressionApplication) = write(io, "\$", _str(comb, latex=true), "\$")


global _lie_expression_index = Dict{LieExpression,Int}()
global _lie_expression_register = Dict{ASCIIString,Tuple{LieExpression,Int}}()

function _get_register_key(ex::LieLinearCombination)
    string('L', join([join([_str_from_objref(x), ':', c]) 
        for (x,c) in sort(ex.terms, 
            lt = (a,b) -> pointer_from_objref(a[1])<pointer_from_objref(b[1])) ],'|'))
end    

function _get_register_key(ex::LieProduct)
    string('P', join([ _str_from_objref(x) for x in ex.factors], ':'))
end

function _get_register_key(ex::LieCommutator)
    string('C', _str_from_objref(ex.A), ':', _str_from_objref(ex.B)) 
end    

function _get_register_key(ex::LieDerivative)
    string('D', _str_from_objref(ex.F)) 
end    

function _get_register_key(ex::LieExponential)
    string('E', _str_from_objref(ex.t),':',_str_from_objref(ex.DF),
           ex.label!=default?string(":",_str_from_objref(ex.label)):"")
end    

function _register(ex::LieExpression)
    key = _get_register_key(ex)
    (ex,i) = get!(_lie_expression_register, key) do
        i = length(_lie_expression_index)
        _lie_expression_index[ex] = i 
        (ex, i)
    end
    ex
end

_str_flat_arg_name(ex::LieExpression) = string("#", _lie_expression_index[ex])

function print_lie_expression_register()
    for (ex, i) in sort(collect(values(_lie_expression_register)), lt = (a,b)-> a[2]<b[2])
        println("#", i, "\t" , _str(ex, flat=true))
     end
end

