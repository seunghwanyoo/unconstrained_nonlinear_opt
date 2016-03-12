function [x] = opt_bfgs(obj,x0,opt,y,h)
% BFGS method for optimization
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

fprintf(' - Running BFGS Method\n');

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

% param for BFGS
fskip = 0;
eps = 10^-8;
Hk = eye(size(hessk));

tic;
k = 0;
if vis > 0
    fprintf('%6s %9s %9s %9s %s\n','iter','f','||grad||','alpha','skip');
    fprintf('%6i %9.2e %9.2e %9.2e %i\n',k,objk,norm(gradk),alpha,fskip);
end
for k = 1:maxiter
    % initialization for alpha
    alpha = 1;
    alpha_u = inf;
    alpha_l = 0;

    %% BFGS Method
    % choose direction
    pk = - Hk*gradk;
    % choose alpha - Wolfe line search
    if ls == 1
        alpha = stepsize_wolfe(alpha,alpha_l,alpha_u,obj,xk,pk,y,h,opt); 
    elseif ls == 2
        alpha = stepsize_backtracking(alpha,obj,xk,pk,y,h,opt);
    end
    % update values
    xk_old = xk;
    xk = xk + alpha*pk;
    gradk_old = gradk;
    gradk = obj.grad(xk,y,h,opt);

    sk = xk - xk_old;
    yk = gradk - gradk_old;
    if k == -1
        Hk = yk'*sk / (yk'*yk) * Hk;
    end
    if sk'*yk > eps*sqrt(sk'*sk)*sqrt(yk'*yk)
        rhok = 1/(yk'*sk);
        Hk = (eye(size(Hk))-rhok*(sk*yk')) * Hk * (eye(size(Hk))-rhok*(yk*sk')) + rhok*(sk*sk');
        fskip = 0;
    else
        fskip = 1;
    end
    if sqrt(gradk'*gradk) < tol %|| abs(objk - objk_old) < 10^(-10)
        fprintf('  converged at %i iter\n',k);
        if vis > 0
            fprintf('%6i %9.2e %9.2e %9.2e %i\n',k,obj.func(xk,y,h,opt),norm(gradk),alpha,fskip);
        end
        break;
    end
        
    if vis > 0
        fprintf('%6i %9.2e %9.2e %9.2e %i\n',k,obj.func(xk,y,h,opt),norm(gradk),alpha,fskip);
    end
end
toc;

% Make sure you set the correct return values
x = xk;
