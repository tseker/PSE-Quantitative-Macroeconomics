We use log-linearization in this homework to solve for the caconical consumer dynamic programming problem, with CRRA utility and with elastic and inelastic labor supply respectively.

`Li_Liang_Seker_loglinear_endogenous_labor.m` is the code for solving the endogenous labor supply case using log linearization. 

`Li_Liang_Seker_loglinear_inelastic_labor.m` is the code for the endogenous labor supply case. We also include a solution using multiple shooting. The result is stark: log linearization generates a very small error.

`solab.m` is the recommended function, attributed to Pual Klain, for solving the dynamic system $Ax_{t+1}=Bx_{t}$ where $A$ is singular.

`cT.m` is the function that calculates the terminal value of consumption to match the steady-state consumption in the multiple shooting algorithm.


