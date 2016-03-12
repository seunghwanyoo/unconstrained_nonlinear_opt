function alpha = stepsize_backtracking(alpha,obj,xk,pk,y,h,opt)
% Line search with Wolfe condition
%
% Input:
%   alpha: current alpha
%   obj: object for function handles
%   xk: current x
%   pk: direction to search
%   y: observed data
%   h: operation matrix
%   opt.c: constraint matrix
%   opt.rho: decreasing factor
% Output:
%   alpha: new alpha - step size
% 
% Author: Seunghwan Yoo

rho = opt.rho;
c1 = 10^(-4);

while obj.func(xk + alpha*pk,y,h,opt) > obj.func(xk,y,h,opt) + c1*alpha*obj.grad(xk,y,h,opt)'*pk
    alpha = rho * alpha;
end