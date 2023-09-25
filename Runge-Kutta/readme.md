
# ODE Solver Library using Runge-Kutta Method (Order 3 and 4)
This is a C library for solving ordinary differential equations (ODEs) of order 3 and 4 using the Runge-Kutta method.

The library provides a set of functions to read system parameters from a file, perform numerical integration, and check convergence conditions.

# Table of Contents
#### Usage
#### Installation
#### File Descriptions

## Usage
To use this library, you can follow these steps:

Clone the repository to your local machine:

Compile the code using a C compiler (e.g., GCC):

gcc seminarski.c rad.c -o ode_solver

### Run the program with the following command:
./ode_solver input.txt p

Replace input.txt with the name of your input data file and p with the order of the Runge-Kutta method (3 or 4).

## Installation
No special installation is required. Simply clone the repository and compile the code as mentioned in the "Usage" section.

## File Descriptions
### seminarski.c: 
The main program that reads input data from a file, performs numerical integration using the Runge-Kutta method of the specified order, and prints the results.

### rad.h: 
Header file containing function declarations and the inclusion guard.

### rad.c: 
Contains the implementation of functions for reading matrices and vectors from a file (ucitaj_matricu_iz_datoteke and ucitaj_vektor_iz_datoteke), performing numerical integration (racunaj), and checking convergence conditions (uslov).


