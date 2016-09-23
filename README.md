# Flows.jl
A Julia package for the symbolic manipulation of flows of nonlinear evolution equations.

Flows.jl was presented at the [CASC 2016](http://www.casc.cs.uni-bonn.de/2016/) workshop
(=>[Slides of the Talk](http://www.harald-hofstaetter.at/Math/SymbolicManipulationOfFlows.pdf))
and is descibed in the paper

>[W. Auzinger](http://www.asc.tuwien.ac.at/~winfried), [H. Hofstätter](http://www.harald-hofstaetter.at), [O. Koch](http://othmar-koch.org), [ Symbolic Manipulation of Flows of Nonlinear Evolution Equations, with Application in the Analysis of Split-Step Time Integrators](http://arxiv.org/pdf/1605.00453.pdf), to appear in [Proceedings of CASC 2016](http://www.casc.cs.uni-bonn.de/2016/).

Flows.jl (or a predecessor written in Perl) was used to carry out and verify the symbolic manipulations 
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
In a Julia notebook type
```julia
Pkg.clone("https://github.com/HaraldHofstaetter/Flows.jl")
```
Note that there is a [new experimental version](https://github.com/HaraldHofstaetter/Flows.jl/tree/lie_experimental), which implements symbolic manipulations of Lie derivatives.
##Examples
To get easy access to the examples, copy them into the home directory:
```julia
cp(joinpath(homedir(), ".julia/v0.4/Flows/examples/"), joinpath(homedir(), "Flows_examples"), remove_destination=true)
```
Then 'Flows_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [Flows.ipynb](https://github.com/HaraldHofstaetter/Flows.jl/blob/master/examples/Flows.ipynb)
+ [DefectLieTrotter.ipynb](https://github.com/HaraldHofstaetter/Flows.jl/blob/master/examples/DefectLieTrotter.ipynb)
+ [ElementaryDifferentials.ipynb](https://github.com/HaraldHofstaetter/Flows.jl/blob/master/examples/ElementaryDifferentials.ipynb)
  
