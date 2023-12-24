We use the value function iteration(VFI) method to solve for the deterministic dynamic programming problem in this project. The main objective is to derive the policy function. We use different flavors of VFI, including non-linear grid points, VFI with linear interpolation, and VFI with golden ratio search. We analyze where the errors for the value function iteration stems from by plotting the Euler Equation Errors against the policy function curvature.

`NGCM_DP_start.m` includes most of the codes for implementation. It includes all algorithims. But figures in question 2 are plotted in `NGCM_DP_Question2.m`.

`NGCM_DP_Question2.m` includes the code for the figures in question 2.

`plotting.m` includes the code for the figures in question 3.

`ValueFunctionIteration_Discrete.m` contains the function for value function iteration with discrete maximizaition. With slight modifications, it can be applied to other scenarios.
