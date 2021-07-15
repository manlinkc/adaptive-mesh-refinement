# Adaptive mesh refinement for the compressible Euler equations using AMReX

Final code for submission for AMReX - written assignment

### AMReX files

**AmrLevelAdv.cpp** – This file is originally from Dr Millmore, it contained C++ code for the flux calculations for solving the advection equation in AMReX. I have modified the code to compute the slope limited MUSCL-Hancock HLLC flux for solving the Euler equations in AMReX. Code which has been modified is identified in the file with comments and my CRSid mc2246.

**Tagging_nd_test2.f90** – This is file is originally from Dr Millmore, it contains Fortran code used to tag cells for refinement. For refining Toro’s Test 2 I have modified the function for tagging cells according to regions of high phi gradient below the threshold. Code which has been modified is identified in the file with comments and my CRSid mc2246.

### Exact Solutions for Toro’s 5 tests

To compute the exact solutions for Toro’s five tests I have created a program which takes in star values for pressure, velocity and density from Table 4.2 in Toro’s textbook ‘Riemann Solvers and Numerical Methods for Fluid Dynamics’ and runs a sampling function to produce to the exact solution.  

**samplingmain.C** – This file reads in inputs from a setting file and reads in x-values from a data file and then samples the exact solution at each of these x-values. The program outputs a data file which contains the x-values and the exact solution at these values. 

**samplingfunctions.C** – This file contains the sampling function and function for reading in parameters from the settings file.

**samplingfunctions.h** – This file contains function headers for the file samplingfunctions.C .

**inputs.txt** – Textfile containing the settings parameters for computing the exact solution for one of Toro’s five tests. 

**master.txt** – Textfile containing the settings parameters for each of Toro’s five tests. The settings for a particular test can be copied from the master.txt file to the inputs.txt file when needed. 

**GNUmakefile** – Builds an executable named ‘mysampling’ which can be compiled by typing make in the directory.   

### Convergence Analysis

The following files are used for the convergence analysis between the 1D numerical and exact solutions. 

**convergence.C** – This file reads in two data files, one which is the exact solution and one which is the numerical solution. It then computes the L1 norm between the exact and numerical solution for density and prints out this value to the terminal. 

**convergenceinputs.txt**  - Textfile where user can specify the name of the files containg the exact and numerical solutions. 
