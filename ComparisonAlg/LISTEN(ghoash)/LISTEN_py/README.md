# LISTEN: Linear Structural Equation Model Learning

A python package for learning linear structural equation models (or equivalently Bayesian networks). This implements the polynomial time
algorithm given in the paper \[1\]. 

# Prerequisites
The package has been tested with Python 2.7 (minor version 12), and depends on the following python packages:

1. numpy,

2. [cvxopt](http://cvxopt.org/) (preferably with GLPK). *Using GLPK is significantly faster than the default solver*.


# Quick start
To learn a SEM from the sample data provided one needs to do the following:
```bash
python code/listen.py data/sample_data_subgaussian.csv weights.csv 0.001 --solver glpk
```
where 0.001 is the regularization parameter. The above command generates a file suffixed by the regularization parameter, i.e., ```weights_lambda_0.001000.csv```. Multiple regularization parameters can be specified at the same time separated by spaces --- one output file is generated per regularization parameter. Omitting the solver option defaults to using the default solver (ConeLP) in cvxopt. Additional options to the program can be found by using the following:
```bash
python code/listen.py -h
```


### References
1. Ghoshal, A. & Honorio, J.. (2018). Learning linear structural equation models in polynomial time and sample complexity. Proceedings of the Twenty-First International Conference on Artificial Intelligence and Statistics (AISTATS), in PMLR 84:1466-1475.

