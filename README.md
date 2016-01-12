# opt_line_search
Matlab implementation of line search methods for nonlinear optimization problem

# implemented methods
1. steepest gradient descent <br>
2. Newton's method <br>
3. modified Newton's method <br>
4. BFGS <br>
5. L-BFGS <br>

# demo program
demo_line_search.m

# options
- method: 1~5 (methods in the above)
- func: test function
- line search: 1 (backtracking), 2 (wolfe)
- rho: parameter for backtracking line search
- m: parameter for L-BFGS
- tol: parameter for stopping criteria



