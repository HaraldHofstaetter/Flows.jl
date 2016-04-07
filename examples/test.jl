using Flows

t = TimeVariable("t")
s = TimeVariable("s")
r = TimeVariable("r")

println("coefficient(t,s) = ", coefficient(t,s))
println("coefficient(t,t) = ", coefficient(t,t))
println("substitute(t,s,r) = ", substitute(t,s,r))
println("substitute(t,t,r) = ", substitute(t,t,r))

ex=TimeLinearCombination(r,1,s,2,t,3)

