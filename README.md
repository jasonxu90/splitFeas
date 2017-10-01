splitFeas
=====

Majorization-minimization for solving split feasibility problems


## Description
This package provides majorization-minimiation (MM) algorithms for solving split feasibility problems. These algorithms were originally proposed for constrained generalized linear model regression [1], and are implemented generally to solve split feasibility problems as detailed in [2]. The goal is to find a point $x \in C$ such that $h(x) \in D$; the sets $C$ and $D$ may be intersections of several closed sets, and need not be convex. The mapping $h$ need not be linear.


## Installation
The package can be installed directly from github using the `devtools` package, which can easily be installed using the command `install.packages("devtools")`.
See https://github.com/hadley/devtools for more details.

To install `splitFeas`, run the following:
```r
library(devtools)
install_github("jasonxu90/splitFeas")
```


## References
1. Xu J, Chi EC, Lange K (2017). "Generalized Linear Model Regression under Distance-to-set Penalties," to appear, *Neural Information Processing Systems*.

2. Xu J, Chi EC, Yang M, Lange K (2017). "A Majorization-minimization Algorithm for Split Feasibility Problems," *arXiv preprint*, https://arxiv.org/abs/1612.05614.

