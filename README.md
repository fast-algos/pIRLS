# pIRLS
This contains a fast IRLS code for solving p-norm regression problems for p at least 2. It is an implementation of the algorithm proposed in the paper, "Fast, Provably convergent IRLS Algorithm for p-norm Linear Regression. Deeksha Adil, Richard Peng and Sushant Sachdeva."

# Algorithm
Algorithm pIRLS is an Iteratively Reweighted Least Squares (IRLS) Algorithm that provably converges for all p at least 2. The algorithm converges geometrically and can thus be used to solve problems to a high accuracy. You may refer to the paper for the analysis and proof of convergence guarantees.

# Using the Code
We have included an implementation in Julia as well as Matlab. The main files have the function implementation. We have two
files, one for graph instances and one for random matrix instances, that can be run directly. For more details on 
these instances and the problems we are solving with them, refer to the paper. The functions can be directly used with other 
inputs as well. Refer to the test files to see how to use them.

---
If you found this code useful in your work, please cite:

```
@incollection{APS19,
title = {Fast, Provably convergent IRLS Algorithm for p-norm Linear Regression},
author = {Adil, Deeksha and Peng, Richard and Sachdeva, Sushant},
booktitle = {Advances in Neural Information Processing Systems 32},
editor = {H. Wallach and H. Larochelle and A. Beygelzimer and F. d\textquotesingle Alch\'{e}-Buc and E. Fox and R. Garnett},
pages = {14166--14177},
year = {2019},
publisher = {Curran Associates, Inc.},
url = {http://papers.nips.cc/paper/9565-fast-provably-convergent-irls-algorithm-for-p-norm-linear-regression.pdf}
}
```
