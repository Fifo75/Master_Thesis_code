# Master_Thesis_code

This file represent the code for the master thesis called "Temporal network, centrality measures, and ordinary differential equations".
The code is written in MATLAB. The code can be find in file called: Thesis_code

The prerequisity for optimal run of the code is to have matlab library called "chebfun".

Description of the files in Thesis_code:
temporal_node_evolution: This function calculates receiver centralities of Katz and matrix exponential and plots the summary plots representing temporal node evolution for each of the centralities.
load_FBsoc_dataset: This function loads the FBsoc dataset and creates the adjacency matrices according to the time granularity (month, week, day).
lagrange_cardinal_basis: Calculates the function handles for the Lagrange Cardinal Basis using Barycentric Interpolation for EQUIDISTANT POINTS.
interpoladion_solver: This function computes the matrix version of the Time-Ordered Exponential based on the section 3.5 using all the theory from chapter 3 from the thesis.
init: This file is base to all the experiments from the thesis. The individual experiments are divided to the sections.
genCoeffMatrix_SP_interval: Computes the coefficient matrix of f(t)*Theta(t-s), with a smooth function f and the Heaviside theta function Theta, in a basis of orthonormal Legendre polynomials.
FBsoc.txt: The text file containing all the data of the FBsoc dataset.
exact_sol: This function conputes exact solution of the Time-Order Exponential using matrix expoential in piece wise constant case or ode45 in continuous case

