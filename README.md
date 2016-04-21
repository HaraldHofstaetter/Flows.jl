# Flows.jl
A Julia package for the symbolic manipulation of flows of nonlinear evolution equations.


Flows.jl (or its Perl predecessor) was used to carry out and verify the symbolic manipulations 
needed for the analysis of error estimators for splitting methods, which is described in the paper

> [W. Auzinger](http://www.asc.tuwien.ac.at/~winfried),
> [H. Hofstätter](http://www.harald-hofstaetter.at),
> [O. Koch](http://othmar-koch.org), 
> [M. Thalhammer](http://techmath.uibk.ac.at/mecht/),
> [Defect-based local error estimators for splitting methods, with application to Schr&ouml;dinger equations, Part III. The nonlinear case](http://www.asc.tuwien.ac.at/preprint/2013/asc19x2013.pdf),
> [J. Comput. and Appl. Math. 273 (2015), pp. 182-204](http://dx.doi.org/10.1016/j.cam.2014.06.012).

Flows.jl consists of only ~1000 lines of Julia code and is fully self-contained, 
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
cp(joinpath(homedir(), ".julia/v0.4/Flows/examples/"), joinpath(homedir(), "Flows_examples"), remove_destination=true)
```
Then 'Flows_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [Flows.ipynb](https://github.com/HaraldHofstaetter/Flows.jl/blob/master/examples/Flows.ipynb)
+ [ElementaryDifferentials.ipynb](https://github.com/HaraldHofstaetter/Flows.jl/blob/master/examples/ElementaryDifferentials.ipynb)
  
