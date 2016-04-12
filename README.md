# Flows.jl
A Julia package for the symbolic manipulations of flows of differential equations.

Flows.jl is a reimplementation of [old Perl code](https://github.com/HaraldHofstaetter/Flows.jl/tree/master/perl), 
which is no longer maintained and is almost impossible to understand.
(The ugly peculiarities of the Perl  syntax were the main motivation to reimplement this stuff in Julia...)

Flows.jl (or its Perl predecessor) was used to carry out and verify the symbolic manipulations 
needed for the analysis of error estimators of splitting methods, which is described in  [these slides for a talk](http://www.harald-hofstaetter.at/Math/Flows_of_Differential_Equations.pdf) and in

> [W. Auzinger](http://www.asc.tuwien.ac.at/~winfried),
> [H. Hofstätter](http://www.harald-hofstaetter.at),
> [O. Koch](http://othmar-koch.org), 
> [M. Thalhammer](http://techmath.uibk.ac.at/mecht/),
> [Defect-based local error estimators for splitting methods, with application to Schr&ouml;dinger equations, Part III. The nonlinear case](http://www.asc.tuwien.ac.at/preprint/2013/asc19x2013.pdf),
> [J. Comput. and Appl. Math. 273 (2015), pp. 182-204](http://dx.doi.org/10.1016/j.cam.2014.06.012).

Flows.jl consists of only ~700 lines of Julia code and is fully self-contained, 
of course with the (very relevant!) exception that it depends on the Julia standard library
(but not on  additional Julia packages).

##Installation
```julia
Pkg.clone("https://github.com/HaraldHofstaetter/Flows.jl")
Pkg.build("Flows")
```
##Examples
To get easy access to the examples, copy them into the home directory:
```julia
cp(joinpath(homedir(), ".julia/v0.4/Flows/examples/"), joinpath(homedir(), "Flows_examples"))
```
Then 'Flows_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [Flows.ipynb](https://github.com/HaraldHofstaetter/Flows.jl/blob/master/examples/Flows.ipynb)
  
