# Unconstrained nonlinear optimization
This is Matlab implementation of unconstrained nonlinear optimization algorithms. The well known six algorithms are implemented. As an application, the constrained least squares (CLS) method for image deconvolution is used. The objective function for CLS is a quadratic function as described below (image deconvolution section). Since it's twice differentiable, we can use Newton's method as well as gradient descent method and quasi-Newton methods.

# Implemented methods
1. Steepest gradient descent <br>
2. Newton's method <br>
3. Modified Newton's method <br>
4. BFGS (Quasi-Newton) <br>
5. L-BFGS (Quasi-Newton) <br>
6. Conjugate gradient <br>

# Image deconvolution
- Degradation model: y = Hx + n
- Deconvolution method: Constrained Least Squares (CLS): <br />
      objective function: min_x 0.5||y-Hx||2 + 0.5*lambda*||Cx||2 <br />

# Description of files
- demo_unconstrained_opt.m: test script
- /opt: includes functions for optimization methods
- /funcs: includes the objective function and its gradient and hessian functions

# Options
- method: 1~6 (methods in the above)
- func: test function
- linesearch: backtracking (1), wolfe (2)
- rho: decreasing factor for backtracking line search
- m: parameter for L-BFGS, number of storage
- tol: parameter for stopping criteria
- maxiter: max iteration
- lambda: param for CLS, regularizing param
- c: param for CLS, regularizing matrix
- vis: param for display, 0:nothing, 1:log, 2:log+figure

# Contact
Seunghwan Yoo (seunghwanyoo2013@u.northwestern.edu)
