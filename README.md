# stokes-dg-figures
Visualisation of simulation data obtained with the ParMooN finite element package in the 
[stokes-dg-experiments project](https://github.com/cristina-v-melnic/stokes-dg-experiments).

## Functionality
- *get_cond_nr.m* -  Computes the system matrix condition number of the linear system of equations corresponding to the Stokes finite-element problem.
- *main.py* - Generates plots for numerical analyses and saves them in specified directories.

By running in *main.py* the functions below, one can obtain the plots in the appendix of the 
[master thesis](https://github.com/cristina-v-melnic/stokes-dg-experiments/blob/main/Master_Thesis_signed.pdf).
- *get_conv_history()* - Collect values of numerical solution errrors from the ParMooN output files from  and plots
the error for every element refinement size, i.e., the convergence history. This is done for 2 examples and multiple
velocity spaces per graph for comparison of performance between numerical schemes. 
- *get_robustness_exp()* - Gets the data for the robustness experiments for the specified example in the directory path
and returns plots of convergence histories for different values of  $Re$ (Reynolds number).
- *get_table()* - Generates text files in a table format containing values for the 'keys' parameters for computations
of different sigmas (parameter of the Stokes DG discretisation) to determine which value would be optimal.

### Examples
#### Norm of error convergence histories
![image](https://user-images.githubusercontent.com/103945852/194314493-8539482b-49f3-4f56-a7ed-d8874ee07370.png)
<p align="center">
 <img src = "https://user-images.githubusercontent.com/103945852/194314201-df4186c0-a79f-4ad3-931e-0a7152575783.png"/>
</p>

#### Pressure-robustness results
![image](https://user-images.githubusercontent.com/103945852/194315365-ba46b034-5080-4ec8-93bd-16cf56d58d6a.png)

#### Parameter grid-search
![image](https://user-images.githubusercontent.com/103945852/194315720-4e28b277-d2fe-4308-9315-3acf5a1bdb86.png)

The figures above are from the contents the 
[master thesis](https://github.com/cristina-v-melnic/stokes-dg-experiments/blob/main/Master_Thesis_signed.pdf), where full description
and discussion of the results can be found.


## Credits
This work was part of my master's thesis project with Prof. Dr. Volker John, Derk Frerichs and Dr. Ulrich Wilbrandt from the
"Numerical Mathematics and Scientific Computing" group at the Weierstrass Institute Berlin. It was financially supported by
the DAAD funding of my master's studies, and the "WIAS Female Master Students Program" funding of my research assistanship.
