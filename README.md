gaussAlg - Algebra of functions for R
=========

[![Build Status](https://travis-ci.org/llaniewski/gaussAlgebra.svg?branch=master)](https://travis-ci.org/llaniewski/gaussAlgebra)

Algebra of functions which looks like:

polynomial * exp( -x^2 )

in high dimensions. And lot of other neat stuff.

Example:
```r
x = seq(-2,2,len=300) # some points for the plot
f = Gauss(0.5) # Gaussian bell with standard deviation 0.5
g = lag( Gauss(0.3), 1) # another Gaussian bell moved to origin 1 (mean=1)
matplot(x, cbind(
  calc(f      ,x), # ploting the values of the function f ...
  calc(g      ,x), # ... and g
  calc(f * g  ,x), # their product, ...
  calc(f + g  ,x), # ... sum ...
  calc(f %% g ,x)  # ... and convolution
), type="l", lty=1, xlab="x", ylab="y", main="Different operations")
legend("topleft",c("f","g","fg (multiplicat)","f + g", "f * g (convolution)"),lty=1, col=1:5)
```
