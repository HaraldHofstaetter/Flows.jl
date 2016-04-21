abstract LieExpression

immutable LieDerivative <: LieExpression
    F::VectorFieldExpression
end

D(F::VectorFieldExpression) = LieDerivative(F)

function _str(DF::LieDerivative; flat::Bool=false, latex::Bool=false) 
    if latex
        string("D_{", _str(DF.F, latex=true), "}")
    else   
        string("D(", _str(DF.F, flat= flat, latex=false), ")")
    end
end    


immutable LieExponential <: LieExpression
    t ::TimeExpression
    DF::LieDerivative
end

exp(t::TimeExpression, DF::LieDerivative) = LieExponential(t, DF)

function _str(E::LieExponential; flat::Bool=false, latex::Bool=false) 
    if latex
        string("\\mathrm{e}^{",  
               typeof(E.t)==TimeLinearCombination?"(":"", 
              _str(E.t, latex=true),
               typeof(E.t)==TimeLinearCombination?")":"", 
              _str(E.DF, latex=true), "}")
    else   
        string("exp(", _str(E.t, flat=flat, latex=false), ",",
               _str(E.DF, flat=flat, latex=false), ")")
    end
end    

immutable LieCommutator <: LieExpression
    A::LieExpression
    B::LieExpression
#    TODO: reactivate constructor once _register is implemented
#    function LieCommutator (A::LieExpression, B::LieExpression)
#        if A==B || A==lie_zero || B==lie_zero
#            lie_zero
#        else
#            _register(new(A,B))
#        end    
#    end        
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
    return LieLinearCombination(Tuple{LieExpression, Real}[(key,val) for (key, val) in d], 0)
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
##immedeately, which is not intended and thus commented out.
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
end

function _str(M::LieProduct; flat::Bool=false, latex::Bool=false) 
    join([string(
        typeof(x)==LieLinearCombination?"(":"", 
        flat?_str_flat_arg_name(x):_str(x, latex=latex), 
        typeof(x)==LieLinearCombination?")":"")
    for x in M.factors], latex?"":"*")
end   

LieExprButNotProduct = Union{LieDerivative,LieExponential,LieLinearCombination, LieCommutator}
#With this trick the construction of Products of Products using only the operator * 
#(not the default constructor) is not possible (i.e. products are flattened).

*(a::LieExprButNotProduct, b::LieExprButNotProduct) = 
    LieProduct(LieExprButNotProduct[a, b])
*(a::LieProduct, b::LieProduct) = LieProduct(vcat(a.factors, b.factors))
*(a::LieExprButNotProduct, b::LieProduct) =
    LieProduct(vcat(a, b.factors))
*(a::LieProduct, b::LieExprButNotProduct) = 
    LieProduct(vcat(a.factors, b))

^(a::LieDerivative, p::Integer)= LieProduct([a for i=1:p])    


string(ex::LieExpression) = _str(ex)
show(io::IO, ex::LieExpression) = print(io, _str(ex))

writemime(io::IO, ::MIME"application/x-latex", ex::LieExpression) = write(io, "\$", _str(ex, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  ex::LieExpression) = write(io, "\$", _str(ex, latex=true), "\$")




immutable LieExSpaceExVarCombination <: SpaceExpression
    lie_ex::LieExpression
    ex::SpaceExpression
    u::SpaceVariable
end
 
combine(lie_ex::LieExpression, ex::SpaceExpression, u::SpaceVariable) =
   LieExSpaceExVarCombination(lie_ex::LieExpression, ex::SpaceExpression, u::SpaceVariable)

combine(lie_ex::LieExpression, u::SpaceVariable) = combine(lie_ex, u, u)

function  combine(lie_ex::LieExpression, Fu::AutonomousFunctionExpression)
   @assert isa(Fu.x, SpaceVariable) string(Fu.fun, "(SpaceVariable) expected")
   combine(lie_ex, Fu, Fu.x)
end

*(lie_ex::LieExpression, u::SpaceVariable) = combine(lie_ex, u)
*(lie_ex::LieExpression, Fu::AutonomousFunctionExpression) = combine(lie_ex, Fu)


function _str(comb::LieExSpaceExVarCombination; flat::Bool=false, latex::Bool=false, arg::SpaceVariable=_no_x_var) 
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
        s = string("combine(",
            _str(comb.lie_ex, flat=flat, latex=false),
            ",",
            _str(comb.ex, flat=flat, latex=false),
            ",",
            _str(comb.u, flat=flat, latex=false),
            ")")
    end
    s
end

string(comb::LieExSpaceExVarCombination) = _str(comb)
show(io::IO, comb::LieExSpaceExVarCombination) = print(io, _str(comb))
writemime(io::IO, ::MIME"application/x-latex", comb::LieExSpaceExVarCombination) = write(io, "\$", _str(comb, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  comb::LieExSpaceExVarCombination) = write(io, "\$", _str(comb, latex=true), "\$")


