# stokes-dg-figures
Visualisation of simulation data obtained with the ParMooN finite element package in the 
[stokes-dg-experiments project](https://github.com/cristina-v-melnic/stokes-dg-experiments).

# Functionality
- *get_cond_nr.m* -  Computes the system matrix condition number of the linear system of equations corresponding to the Stokes finite-element problem.
- *main.py* - Generates plots for numerical analyses and saves them in specified directories.

By running in *main.py* the functions below, one can obtain the plots in the appendix of the 
[master thesis](https://github.com/cristina-v-melnic/stokes-dg-experiments/blob/main/Master_Thesis_signed.pdf).
- *get_conv_history()* - Collect values of numerical solution errrors from the ParMooN output files from  and plots
the error for every element refinement size, i.e., the convergence history. This is done for 2 examples and multiple
velocity spaces per graph for comparison of performance between numerical schemes. 
- *get_robustness_exp()* - Gets the data for the robustness experiments for the specified example in the directory path
and returns plots of convergence histories for different values of  $Re$ (Reynolds number).
- *get_table()* - Generates text files in a latex table format containing values for the 'keys' parameters for computations
of different sigmas (parameter of the Stokes DG discretisation) to determine which value would be optimal.


##### Credits
This work was part of my master's thesis project with Prof. Dr. Volker John, Derk Frerichs and Dr. Ulrich Wilbrandt from the
"Numerical Mathematics and Scientific Computing" group at the Weierstrass Institute Berlin. It was financially supported by
the DAAD funding of my master's studies, and the "WIAS Female Master Students Program" funding of my research assistanship.
