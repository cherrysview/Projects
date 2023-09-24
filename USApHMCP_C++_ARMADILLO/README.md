# GVNS (General Variable Neighborhood Search) for UFLP (Uncapacitated Facility Location Problem)

This C++ code implements the General Variable Neighborhood Search (GVNS) algorithm to solve the Uncapacitated Facility Location Problem (UFLP). The code uses the Armadillo library for linear algebra operations.

## Introduction
The Uncapacitated Facility Location Problem (UFLP) is a classic combinatorial optimization problem where you need to decide which facilities to open among a set of potential facility locations to minimize the total cost, including both facility opening costs and assignment costs for customers to open facilities.

## Prerequisites
Before running the code, make sure you have the Armadillo library installed. You can download Armadillo from http://arma.sourceforge.net/ and follow the installation instructions.

## Input Data
The input data for the UFLP problem is read from a file named ulaz_UFLP.txt. The format of the input file is as follows:

scss
Copy code
n (number of nodes)
p (number of facilities to be opened)
d11,d12,...,d1n (distance matrix)
d21,d22,...,d2n
...
dn1,dn2,...,dnn
c11,c12,...,c1n (cost matrix)
c21,c22,...,c2n
...
cn1,cn2,...,cnn
alpha (discount factor)
beta (the deadline travelling time)
opt_sol (known optimal solution, if available or the best feasible solution from literature)

## Algorithm Description
The code implements the GVNS algorithm, which is a metaheuristic for solving combinatorial optimization problems. It consists of the following main steps:

 Initialize: Generate an initial random solution (set of open facilities).

 Shake: Perturb the current solution by randomly closing some facilities and opening new ones.

Local Search (VND): Perform a variable neighborhood descent (VND) on the shaken solution to improve it. Improvement based on demand
for each non-hub node reallocated to other hub (fixed hub location) evaluate improvement for total demand.

Update: If the improved solution is better than the current best solution, update the best solution.
Repeat: Repeat steps 2-4 for a predefined number of iterations.

The algorithm aims to find a good solution to the USApHMCP by iteratively exploring different neighborhoods of solutions.

## Running the Code
Compile the code with your C++ compiler, ensuring that you have the Armadillo library properly linked. You can use the following command:

The program will read the input data from the ulaz_USApHMCP.txt file and output various results, including:

Average initial solution time (t)
Average total solution time (t_tot)
Best solution found
Average gap between the solutions found and the known optimal solution (if available)
Standard deviation of the gap

## Notes
The code uses randomization to generate initial solutions and explore neighborhoods, so it may produce different results in each run.

The performance of the GVNS algorithm depends on various factors, including the problem instance and parameter settings. 

You may need to adjust parameters like maxiter, kmax, and lmax for different problem sizes.

The code calculates and reports the gap between the solutions found and the known optimal solution (if available). 
A smaller gap indicates a better-quality solution.
