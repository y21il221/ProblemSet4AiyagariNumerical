1. "aiyagari_origin.m", iteration without using any technique, K*=30.278, r*=1.0101, Gini=0.315, runtime is 363.97 seconds.

2. "aiyagari_policy_iteration.m", iteration with policy iteraction, K*=30.278, r*=1.0101, Gini=0.315, runtime is 359.95 seconds.

3. "aiyagari_inter_linear.m", for the linear interpolation, first use the vectorized maximization to find out 80 grid points. Then use the linear interpolation with 800 hundreds grid points, and use vectorized maximiztion to find out the maximizer a' for each of the 800 values of a. K*=30.691, r*=1.0098, Gini=0.2447, runtime is 159.89 seconds.


4. "aiyagari_inter_cubic.m", for the cubic interpolation, first use the vectorized maximization to find out 80 grid points. Then use the cubic interpolation with 80 hundreds grid points, and use the GoldenMethodSearch to find out the maximizer a' for each of the 80 values of a for the continuous function, this gives another 80 new points for the cubic interpolation, iteract until the new value function is close to the previous one. Then plug 800 hundreds grid points into the cubic interpolation, and use vectorized maximiztion to find out the maximizer a' for each of the 800 values of a. K*=30.687, r*=1.0098, Gini=0.2472, runtime is 288.86 seconds

5. Compared to the Huggett model, the wealth distribution of the Aiyagari model has long right tail.

6. The graphs of policy functions, wealth distribution, and the Lorenz curve of the cubic interpolation are attached.

7. Please use the "GoldenSectionSearch.m" uploaded, some editions have been made.