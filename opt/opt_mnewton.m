function [x] = opt_mnewton(obj,x0,opt,y,h)
% Modified Newton's method for optimization
% Inputs:
%   obj: functions to evaluate objective and its gradient and hessian
%   x0:  initial point
%   opt: algorithmic options
%   y: observed data
%   h: operation matrix
%   c: constraint matrix
% Output:
%   sol: solution - converged point
%
% Author: Seunghwan Yoo

fprintf(' - Running Modified Newton''s Method\n');

% k-th (k=0) function, gradient, hessian
objk  = obj.func(x0,y,h,opt);
gradk = obj.grad(x0,y,h,opt);
hessk = obj.hess(x0,y,h,opt);
xk = x0;

% options
disp(opt);
ls = opt.linesearch;
tol = opt.tol;
maxiter = opt.maxiter;
vis = opt.vis;

% param for line search
alpha = 1; % initial value
inf = 10^10; % value for infinity

tic;
k = 0;
if vis > 0
    fprintf('%6s %9s %9s %9s %s\n','iter','f','||grad||','alpha');
    fprintf('%6i %9.2e %9.2e %9.2e %i\n',k, objk, norm(gradk), alpha);
end
for k = 1:maxiter
    % initialization for alpha
    alpha = 1;
    alpha_u = inf;
    alpha_l = 0;

    %% Modified Newton's Method    
    % Cholesky with added multiple of the identity
    beta = 10^(-3);
    d_min = min(diag(hessk));
    if d_min > 0
        tao = 0;
    else
        tao = -d_min + beta;
    end
    while 1
        RTR = hessk + tao*eye(size(hessk));
        [R,p] = chol(RTR);
        if p == 0
            break;
        else
            tao = max(2*tao,beta);
        end
    end
    %Bk = R'*R;  % Bk * pk = -gradk
    % choose direction
    pk = - R \ (R' \ gradk); % instead of  pk = -Bk \ gradk;
    % choose step size
    if ls == 1
        alpha = stepsize_wolfe(alpha,alpha_l,alpha_u,obj,xk,pk,y,h,opt);
    elseif ls == 2
        alpha = stepsize_backtracking(alpha,obj,xk,pk,y,h,opt);
    end
    % update values
    xk_old = xk;
    xk = xk + alpha*pk;
    if norm(xk-xk_old)/norm(xk_old) < tol
        fprintf('  converged at %i iter\n',k);
        if vis > 0
            fprintf('%6i %9.2e %9.2e %9.2e %i\n',k,obj.func(xk,y,h,opt),norm(gradk),alpha);
        end
        break;
    end
    gradk = obj.grad(xk,y,h,opt);
    hessk = obj.hess(xk,y,h,opt);
        
    if vis > 0
        fprintf('%6i %9.2e %9.2e %9.2e %i\n',k,obj.func(xk,y,h,opt),norm(gradk),alpha);
    end
end
toc;

% Make sure you set the correct return values
x = xk;
