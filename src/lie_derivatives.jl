Operator = AutonomousFunction

immutable OperatorLinearCombination 
    terms :: Array{Tuple{Union{Operator,OperatorLinearCombination}, Real},1}
    function OperatorLinearCombination(terms
        :: Array{Tuple{Union{Operator,OperatorLinearCombination}, Real},1}, dummy::Int)
       # dummy only to make it distinguishable from the constructor below
       new(terms)
   end    
end

OperatorExpression = Union{Operator,OperatorLinearCombination}

function OperatorLinearCombination(terms :: Array{Tuple{OperatorExpression, Real},1})
    d = Dict{OperatorExpression,Real}()
    for (x,c) in terms
        #@assert isa(c, Real) "Real expected"
        if isa(x, OperatorLinearCombination)
            # Nested OperatorLinearCombination are expanded into parent OperatorLinearCombination
            for (x1, c1) in x.terms
                get!(d, x1, 0) 
                d[x1] += c*c1
            end
        else
            #@assert isa(x, OperatorExpression) "OperatorExpression expected"
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
    return OperatorLinearCombination(Tuple{OperatorExpression, Real}[(key,val) for (key, val) in d], 0)
end

OperatorLinearCombination(x...) = OperatorLinearCombination(Tuple{OperatorExpression, Real}[ (x[i],x[i+1]) for i=1:2:length(x) ])

+(a::OperatorExpression, b::OperatorExpression) = OperatorLinearCombination(a,1, b, 1)
-(a::OperatorExpression, b::OperatorExpression) = OperatorLinearCombination(a,1, b,-1)
-(a::OperatorExpression) = (-1)*a
+(a::OperatorExpression) = a
*(f::Real, ex::OperatorExpression) = OperatorLinearCombination(ex,f)
*(f::Real, ex::OperatorLinearCombination) = OperatorLinearCombination( Tuple{OperatorExpression, Real}[ (x, f*c) for (x, c) in ex.terms ] )
*(ex::OperatorExpression, f::Real) = f*ex


function _str(ex::OperatorLinearCombination; flat::Bool=false, latex::Bool=false) 
    if length(ex.terms) == 0 
        return latex?"0":"t_zero"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)==OperatorLinearCombination?"(":"", 
            flat?_str_flat_arg_name(x):_str(x, latex=latex), 
            typeof(x)==OperatorLinearCombination?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   

string(ex::OperatorLinearCombination) = _str(ex)
show(io::IO, ex::OperatorLinearCombination) = print(io, _str(ex))

writemime(io::IO, ::MIME"application/x-latex", ex::OperatorLinearCombination) = write(io, "\$", _str(ex, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  ex::OperatorLinearCombination) = write(io, "\$", _str(ex, latex=true), "\$")


coefficient(ex::Operator, v::Operator) = (ex==v ? 1 : 0)

function coefficient(ex::OperatorLinearCombination, v::Operator)
    c = 0
    for (ex1, c1) in ex.terms
       c += c1*coefficient(ex1, v)
    end
    c
end 

abstract LieExpression

immutable LieDerivative <: LieExpression
    F::OperatorExpression
end

D(F::OperatorExpression) = LieDerivative(F)

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


immutable LieMonomial <: LieExpression
    factors::Array{Union{LieDerivative,LieExponential},1}
end

*(a::Union{LieDerivative,LieExponential}, b::Union{LieDerivative,LieExponential}) = 
    LieMonomial(Union{LieDerivative,LieExponential}[a, b])
*(a::LieMonomial, b::LieMonomial) = LieMonomial(vcat(a.factors, b.factors))
*(a::Union{LieDerivative,LieExponential}, b::LieMonomial) =
    LieMonomial(vcat(a, b.factors))
*(a::LieMonomial, b::Union{LieDerivative,LieExponential}) = 
    LieMonomial(vcat(a.factors, b))

function _str(M::LieMonomial; flat::Bool=false, latex::Bool=false) 
    join([_str(x, flat=flat, latex=latex) for x in M.factors], latex?"":"*")
end    


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

function *(a::LieLinearCombination, b::LieLinearCombination)
    LieLinearCombination( 
        reshape(Tuple{LieExpression, Real}[ (x*y, c*d)   
                for (x,c) in a.terms, (y,d) in b.terms], 
                length(a.terms)*length(b.terms)) )
end   

function *(a::Union{LieDerivative,LieExponential,LieMonomial}, b::LieLinearCombination)
    LieLinearCombination( Tuple{LieExpression, Real}[ (a*x, c)  for (x,c) in b.terms] ) 
end   
function *(a::LieLinearCombination, b::Union{LieDerivative,LieExponential,LieMonomial})
    LieLinearCombination( Tuple{LieExpression, Real}[ (x*b, c)  for (x,c) in a.terms] ) 
end   


function _str(ex::LieLinearCombination; flat::Bool=false, latex::Bool=false) 
    if length(ex.terms) == 0 
        return latex?"0":"t_zero"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)==LieLinearCombination?"(":"", 
            flat?_str_flat_arg_name(x):_str(x, latex=latex), 
            typeof(x)==LieLinearCombination?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   



string(ex::LieExpression) = _str(ex)
show(io::IO, ex::LieExpression) = print(io, _str(ex))

writemime(io::IO, ::MIME"application/x-latex", ex::LieExpression) = write(io, "\$", _str(ex, latex=true), "\$")
writemime(io::IO, ::MIME"text/latex",  ex::LieExpression) = write(io, "\$", _str(ex, latex=true), "\$")


