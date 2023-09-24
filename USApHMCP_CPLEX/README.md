## CPLEX ILOG Project Readme
## Introduction
This project is an implementation of a mathematical optimization problem using the CPLEX ILOG libraries. 
The code solves a specific problem known as the USApHMCP (Uncapacitated Single Allocation p-Hub Maximal Covering Problem) using a integer linear programming (ILP) approach.

## Problem Description
The aim of model is to locate the hub, and to allocate non-hub nodes to the located hub nodes where each node can be potential hub-node. 
The hub can maximize the demand covered by deadline traveling time.

## Getting Started
Before you can run the code, you need to set up the environment and ensure you have the necessary input data.

## Installation: 
Make sure you have the CPLEX ILOG libraries installed on your system. You can obtain these libraries from IBM CPLEX.

## Input Data: 
The code expects input data to be in a file named "USApHMCP_ulaz.txt." Ensure that this file contains the required problem data in the specified format.

## Compilation: 
Compile the code using a C++ compiler. Ensure that you link it with the CPLEX ILOG libraries properly.

## Execution: 
Run the compiled executable. The code will read the input data, formulate the ILP problem, and use the CPLEX solver to find the optimal solution.

## Input Data Format
The input data should be in the "USApHMCP_ulaz.txt" file and should be formatted as follows:

The first two lines contain two integers, n and p, representing problem parameters number of nodes and number of hubs.
The next n lines contain n floating-point values, representing the matrix D.
The next n lines contain n floating-point values, representing the matrix C.
The following lines contain floating-point values for alpha, beta, and integer 777 means check.
Please ensure that your input data is formatted correctly according to the problem requirements.

## Code Structure
The code is structured as follows:

It reads the input data and initializes the necessary variables and matrices.
It formulates the ILP problem using CPLEX ILOG constructs, defining decision variables, constraints, and the objective function.
It sets various CPLEX solver parameters to control the optimization process.
It solves the ILP problem using the CPLEX solver.
It prints the results, including the status, objective function value, execution time, and other relevant information, to an output file named "izlaz_USApHMCP.txt."

## Output
After running the code, you can find the results in the "izlaz_USApHMCP.txt" file. The output includes the solution status, objective function value, execution time, and other useful information about the optimization process.

## Additional Notes
Make sure to customize the input data file name or provide the correct instance data in the specified format.
The code includes comments to explain various parts of the implementation in Serbian, making it easier to understand and modify if needed.
