function grad = func3_grad(x,y,h,opt)
% Gradient of the obj. function
% Input: 
%   x: original image (vector)
%   y: degraded image (vector)
%   h: impulse response (lexicographically arranged)
%   opt: option for optimization - lambda, c
% Output:
%   grad: gradient

lambda = opt.lambda;
c = opt.c;

grad = (h)'*h*x - (h)'*y + lambda*(c)'*c*x;
