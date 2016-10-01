#This magic code stolen from SymPy, cf.
#https://github.com/jverzani/SymPy.jl/blob/master/src/utils.jl

macro t_vars(x...) 
    q=Expr(:block)
    if length(x) == 1 && isa(x[1],Expr)
        @assert x[1].head === :tuple "@t_vars expected a list of symbols"
        x = x[1].args
    end
    for s in x
        @assert isa(s,Symbol) "@t_vars expected a list of symbols"
        push!(q.args, Expr(:(=), s, Expr(:call, :TimeVariable, Expr(:quote, string(s)))))
    end
    push!(q.args, Expr(:tuple, x...))
    eval(Main, q)
end

macro x_vars(x...) 
    q=Expr(:block)
    if length(x) == 1 && isa(x[1],Expr)
        @assert x[1].head === :tuple "@x_vars expected a list of symbols"
        x = x[1].args
    end
    for s in x
        @assert isa(s,Symbol) "@x_vars expected a list of symbols"
        push!(q.args, Expr(:(=), s, Expr(:call, :SpaceVariable, Expr(:quote, string(s)))))
    end
    push!(q.args, Expr(:tuple, x...))
    eval(Main, q)
end

macro funs(x...) 
    q=Expr(:block)
    if length(x) == 1 && isa(x[1],Expr)
        @assert x[1].head === :tuple "@funs expected a list of symbols"
        x = x[1].args
    end
    for s in x
        @assert isa(s,Symbol) "@funs expected a list of symbols"
        #push!(q.args, Expr(:(=), s, Expr(:call, :AutonomousFunction, Expr(:quote, string(s)))))
        push!(q.args, Expr(:(=), s, Expr(:call, :VectorFieldVariable, Expr(:quote, string(s)))))
    end
    push!(q.args, Expr(:tuple, x...))
    eval(Main, q)
end

macro nonautonomous_funs(x...) 
    q=Expr(:block)
    if length(x) == 1 && isa(x[1],Expr)
        @assert x[1].head === :tuple "@nonautonomous_funs expected a list of symbols"
        x = x[1].args
    end
    for s in x
        @assert isa(s,Symbol) "@nonautonomous_funs expected a list of symbols"
        push!(q.args, Expr(:(=), s, Expr(:call, :NonAutonomousFunction, Expr(:quote, string(s)))))
    end
    push!(q.args, Expr(:tuple, x...))
    eval(Main, q)
end


