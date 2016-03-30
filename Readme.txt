The cubic iteration code is rewritten, the code mainly includes 5 parts,

1. Coarse grid search with 51 grids

2. Cubic interpolate the value function by the 51 grids, use the golden section search to maximize the value function, and interpolate the value function by the new 51 values, if the new values are close to the original values of the value function, the loop ends. If the values are not close to the original values, go back to step 2

3. Use the interpolate function to guess the continuous value of the 801 grids, and implement the grid search

4. Cubic interpolate the value function by the 801 grids, use the golden section search to maximize the value function, and interpolate the value function by the new 801 values, if the new values are close to the original values of the value function, the loop ends. If the values are not close to the original values, go back to step 4

5. Use the method you suggest to adjust the policy function and calculate the invariant distribution





K*=30.6866, r*=1.0098, Gini=0.2125, runtime is 1380 seconds

The graphs of the policy function, wealth distribution, and Lorenz curve are attached

The GoldenSectionSearch m-file is modified