function x = opt_cg(obj,x0,opt,y,h)
% Conjugate gradient method for optimization 
%  ***only for quadratic problem, 1/2 x'Ax - bx (=> Ax = b)
% Inputs:
%   obj: functions to evaluate objective and its gradient and hessian
%   x0:  initial point
%   opt: algorithmic options
%   y: observed data
%   h: operation matrix
%   opt.c: constraint matrix
% Output:
%   sol: solution - converged point
%
% Author: Seunghwan Yoo

fprintf(' - Running Conjugate Gradient Method\n');

% k-th (k=0) function, gradient, hessian
%objk  = obj.func(x0,y,h,opt);
gradk = obj.grad(x0,y,h,opt);
hessk = obj.hess(x0,y,h,opt);
x = x0;

%Initialization
tol = opt.tol;

A = hessk;
r = gradk; %A*x-b;
p = -r;
k = 0;
abs_r = sqrt(r'*r);

tic;
while abs_r > tol
    k = k+1;
    alpha = (r'*r)/(p'*A*p);
    x = x + alpha*p;
    r_old = r;
    r = r + alpha*A*p;
    beta = (r'*r)/(r_old'*r_old);
    p = -r + beta*p;
    abs_r = sqrt(r'*r);
end
fprintf('  converged at %i iter\n',k);
toc;
