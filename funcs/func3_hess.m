function H = func3_hess(x,y,h,opt)
% Hessian of the objective function
% Input: 
%   x: original image (vector)
%   y: degraded image (vector)
%   h: impulse response (lexicographically arranged)
%   opt: option for optimization - lambda, c
% Output:
%   H: Hessian

lambda = opt.lambda;
c = opt.c;
hth = h'*h;

H = hth + lambda*(c)'*c;