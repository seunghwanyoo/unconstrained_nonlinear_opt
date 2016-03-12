function [x] = opt_lbfgs(obj,x0,opt,y,h)
% L-BFGS method for optimization
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

fprintf(' - Running L-BFGS Method\n');

% k-th (k=0) function, gradient, hessian
objk  = obj.func(x0,y,h,opt);
gradk = obj.grad(x0,y,h,opt);
%hessk = obj.hess(x0,y,h,opt);
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

% param for L-BFGS
m = opt.m;
n = length(x0);
si = zeros(n,m);
yi = zeros(n,m);
rhoi = zeros(1,m);
ai = zeros(1,m);

tic;
k = 0;
if vis > 0
    fprintf('%6s %9s %9s %9s\n','iter','f','||grad||','alpha');
    fprintf('%6i %9.2e %9.2e %9.2e %i\n',k,objk,norm(gradk),alpha);
end
for k = 1:maxiter
    % initialization for alpha
    alpha = 1;
    alpha_u = inf;
    alpha_l = 0;

    %% L-BFGS Method
    % Hk0
    if k == 1
        Hk0 = eye(n);
    else
        gammak = sk'*yk / (yk'*yk);%sk_old'*yk_old / (yk_old'*yk_old);
        Hk0 = gammak*eye(n);
    end

    % compute pk = - Hk*gradk;
    l = 0;
    q = gradk;
    if k-1 >= 1
        for i = k-1:-1:max(k-m,1) % k-1 to k-m
            idx = m - l;
            if idx == m
                si(:,idx) = sk; %xk - xk_old;
                yi(:,idx) = yk; %gradk - gradk_old;
                rhoi(idx) = 1/(yi(:,idx)'*si(:,idx));
            end
            ai(idx) = rhoi(idx)*si(:,idx)'*q;
            q = q - ai(idx)*yi(:,idx);
            l = l + 1;
        end
    end
    r = Hk0*q;
    if k-1 >= 1
        for i = max(k-m,1):1:k-1
            l = l - 1;
            idx = m - l;
            beta = rhoi(idx)*yi(:,idx)'*r;
            r = r + si(:,idx)*(ai(idx)-beta);
        end
    end
    pk = - r;

    % compute alpha
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

    for i = 1:m-1
        si(:,i) = si(:,i+1);
        yi(:,i) = yi(:,i+1);
        rhoi(i) = rhoi(i+1);
    end

    if sqrt(gradk'*gradk) < tol %|| abs(objk - objk_old) < 10^(-10)
        fprintf('  converged at %i iter\n',k);
        if vis > 0
            fprintf('%6i %9.2e %9.2e %9.2e %i\n',k,obj.func(xk,y,h,opt),norm(gradk),alpha);
        end
        break;
    end
        
    if vis > 0
        fprintf('%6i %9.2e %9.2e %9.2e %i\n',k,obj.func(xk,y,h,opt),norm(gradk),alpha);
    end
end
toc;

% Make sure you set the correct return values
x = xk;
