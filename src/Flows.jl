module Flows

import 
    Base: (*), +, -, string, show

export TimeExpression, TimeVariable, TimeLinearCombination
export coefficient, substitute, simplify


###################################################################################################

abstract TimeExpression

type TimeVariable <: TimeExpression
   name::AbstractString
   latex::AbstractString
end 

function TimeVariable(name::AbstractString; 
                      latex::AbstractString=name)
    TimeVariable(name, latex)
end   

string(t::TimeVariable) = t.name
show(io::IO, t::TimeVariable) = print(io, string(t))

type TimeLinearCombination <: TimeExpression
   terms :: Array{Tuple{TimeExpression, Real},1}
end

TimeLinearCombination(x...) = simplify(TimeLinearCombination([ (x[i],x[i+1]) for i=1:2:length(x) ]))

+(a::TimeExpression, b::TimeExpression) = TimeLinearCombination(a,1, b, 1)
-(a::TimeExpression, b::TimeExpression) = TimeLinearCombination(a,1, b,-1)
-(a::TimeExpression) = (-1)*a
*(f::Real, ex::TimeVariable) = TimeLinearCombination(ex,f)
*(f::Real, ex::TimeLinearCombination) = TimeLinearCombination( [ (x, f*c) for (x, c) in ex.terms ] )
*(ex::TimeExpression, f::Real) = f*ex

function string(ex::TimeLinearCombination) 
    if length(ex.terms) == 0 
        return "0"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)!=TimeVariable?"(":"", 
            string(x), 
            typeof(x)!=TimeVariable?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   

show(io::IO, t::TimeLinearCombination) = print(io, string(t))

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
    TimeLinearCombination([(key,val) for (key, val) in d])
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


###################################################################################################

type AutonomousFunction
   name::AbstractString
   latex::AbstractString
end

function AutonomousFunction(name::AbstractString; 
                            latex::AbstractString=name)
    AutonomousFunction(name, latex)
end   

string(t::TimeVariable) = t.name
show(io::IO, t::TimeVariable) = print(io, string(t))

###################################################################################################
# Essentially the same stuff again with Time substituted by Space...

abstract SpaceExpression

type SpaceTimeVariable <: SpaceExpression
   name::AbstractString
   latex::AbstractString
end 

function SpaceVariable(name::AbstractString; 
                      latex::AbstractString=name)
    SpaceVariable(name, latex)
end   

string(t::SpaceVariable) = t.name
show(io::IO, t::SpaceVariable) = print(io, string(t))


type SpaceLinearCombination <: SpaceExpression
   terms :: Array{Tuple{SpaceExpression, Real},1}
end

SpaceLinearCombination(x...) = simplify(SpaceLinearCombination([ (x[i],x[i+1]) for i=1:2:length(x) ]))

+(a::SpaceExpression, b::SpaceExpression) = SpaceLinearCombination(a,1, b, 1)
-(a::SpaceExpression, b::SpaceExpression) = SpaceLinearCombination(a,1, b,-1)
-(a::SpaceExpression) = (-1)*a
*(f::Real, ex::SpaceVariable) = SpaceLinearCombination(ex,f)
*(f::Real, ex::SpaceLinearCombination) = SpaceLinearCombination( [ (x, f*c) for (x, c) in ex.terms ] )
*(ex::SpaceExpression, f::Real) = f*ex

function string(ex::SpaceLinearCombination) 
    if length(ex.terms) == 0 
        return "0"  #empty linear combination
    else    
        s = join([join([c>=0?"+":"-", abs(c)==1?"":abs(c),
            typeof(x)!=SpaceVariable?"(":"", #TODO: which other types don't need parantheses?
            string(x), 
            typeof(x)!=SpaceVariable?")":"", 
        ]) for (x, c) in ex.terms])
        return s[1]=='+' ? s[2:end] : s
    end    
end   

show(io::IO, t::SpaceTimeLinearCombination) = print(io, string(t))


substitute(ex::SpaceVariable, this::SpaceVariable, by::SpaceExpression) = (ex==this ? by : ex)
substitute(ex::SpaceVariable, this::TimeVariable, by::TimeExpression) = ex 

function substitute(ex::SpaceLinearCombination, this::SpaceVariable, by::SpaceExpression)
    TimeLinearCombination([(substitute(x, this, by), c) for (x, c) in ex.terms])
end

function substitute(ex::SpaceLinearCombination, this::ThisVariable, by::ThisExpression)
    TimeLinearCombination([(substitute(x, this, by), c) for (x, c) in ex.terms])
end



###################################################################################################




end #module Flows
