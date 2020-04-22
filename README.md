# NaturalES.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.org/francescoalemanno/NaturalES.jl.svg?branch=master)](https://travis-ci.com/francescoalemanno/NaturalES.jl)
[![codecov.io](http://codecov.io/github/francescoalemanno/NaturalES.jl/coverage.svg?branch=master)](http://codecov.io/github/francescoalemanno/NaturalES.jl?branch=master)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://francescoalemanno.github.io/NaturalES.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://francescoalemanno.github.io/NaturalES.jl/dev)
-->

This package implements the optimization methods described in
[_Wierstra, et al "Natural Evolution Strategies", JMLR (2014)_](http://www.jmlr.org/papers/volume15/wierstra14a/wierstra14a.pdf).
this implementation follows the KISS™ principle, it can be used as

# Usage

```julia
function rosenbrock(x::AbstractVector{T}) where T
    s=(1.0 - x[1])^2
    for i in 1:(length(x)-1)
        s+=100.0 * (x[i+1] - x[i]^2)^2
    end
    return s
end

optimize(rosenbrock,[0.3,0.6],1.0,sNES) # separable natural es.

(sol = [0.9999902815083116, 0.9999805401026993], cost = 9.450201922031972e-11)


optimize(rosenbrock,[0.3,0.6],1.0,xNES) # exponential natural es.

(sol = [0.9999999934969991, 0.9999999871800216], cost = 4.574949214506023e-17)
```

for further info in Julia type `?optimize`.

# Tips:

* Use xNES for hard problems with strongly correlated variables
* Use sNES for high dimensional problems that exhibit many local minima
* Use sNES for problems with mostly separable variables

## Other packages

look at the excellent `BlackBoxOptim`, or `Optim`
