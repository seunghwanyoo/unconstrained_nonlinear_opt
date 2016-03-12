% This is Matlab implementation of unconstrained nonlinear optimization 
% algorithms. The well known six algorithms are implemented. As an 
% application, the constrained least squares (CLS) method for image
% deconvolution is used. The objective function for CLS is a quadratic 
% function as described below (image deconvolution section). Since it's 
% twice differentiable, we can use Newton's method as well as gradient 
% descent method and quasi-Newton methods.
%
% Implemented methods
% 1. Steepest gradient descent 
% 2. Newton's method 
% 3. Modified Newton's method 
% 4. BFGS (Quasi-Newton) 
% 5. L-BFGS (Quasi-Newton) 
% 6. Conjugate gradient 
%
% Degradation model: y = Hx + n 
% Deconvolution:
%   CLS solver: min_x 0.5||y-Hx||2 + 0.5*lambda*||Cx||2
%
% Author: Seunghwan Yoo (seunghwanyoo2013@u.northwestern.edu)

clear; close all;
addpath(genpath('.'));

% param.opt
% 1. Steepest gradient descent 
% 2. Newton's method 
% 3. Modified Newton's method 
% 4. BFGS (Quasi-Newton) 
% 5. L-BFGS (Quasi-Newton) 
% 6. Conjugate gradient 
param.opt = 6;
param.blur = 1; % 1:Gaussian kernel, 2:User defined

%% original image
x0_whole = im2double(imread('peppers.png'));%'greens.jpg';
if ndims(x0_whole) > 1
    x0_whole = rgb2gray(x0_whole);
end
x_2d = x0_whole(201:220,201:220); % original image (20x20)
x = x_2d(:); % vectorized

%% blur kernel & laplacian kernel
switch (param.blur)
    case 1
        h0_2d = fspecial('gaussian',[11,11],2);
    case 2
        h0_2d = [1 1 1; 1 1 1; 1 0 0];%h0 = ones(5,5); % blur kernel
        h0_2d = h0_2d/sum(sum(h0_2d)); % blur kernel
end
c0_2d = [0 0.25 0; 0.25 -1 0.25; 0 0.25 0]; % 2D Laplacian for CLS

%% create operator matrix for lexicographic notation
tic; [h,h_2d] = create_lexicoH(x_2d,h0_2d); toc;
tic; [c,c_2d] = create_lexicoH(x_2d,c0_2d); toc;

%% degradation
fprintf('\n== Degradation\n');
y_b = h*x;                            % blurred image (vector)
y_2d_b = reshape(y_b,size(x_2d));     % blurred image (2D);
y_2d = imnoise(y_2d_b,'gaussian',0,0.01);  % noisy image (2D);
y = y_2d(:);                               % noisy image (vector);
figure, imshow(x_2d); title('original');
figure, imshow(y_2d); title('degraded (blur+noise)');


%% non-blind deconvolution (with noise, known y,h, get x)
%%% CLS
fprintf('== CLS with optimization methods\n');
opt.linesearch = 1; % 1:wolfe, 2:backtracking
opt.rho = 0.5;      % param for backtracking line search
opt.tol = 10^(-5);  % param for stopping criteria
opt.maxiter = 10^3; % param for max iteration
opt.m = 6;          % param for L-BFGS, num for storage
opt.lambda = 10;    % param for CLS, regularizing param
opt.c = c;          % param for CLS, regularizing matrix
opt.vis = 0;        % param for display, 0:nothing,1:log,2:log+figure
obj.func = @func3;  % func1:LS, func2:CLS w/ I, func3:CLS w/ C
obj.grad = @func3_grad;
obj.hess = @func3_hess;
x0 = y;
switch (param.opt)
    case 1
        [x_cls_i] = opt_gd(obj,x0,opt,y,h);
        figure, imshow(reshape(x_cls_i,size(x_2d))); title('gradient descent');
    case 2
        [x_cls_i] = opt_newton(obj,x0,opt,y,h);
        figure, imshow(reshape(x_cls_i,size(x_2d))); title('Newtons');
    case 3
        [x_cls_i] = opt_mnewton(obj,x0,opt,y,h);
        figure, imshow(reshape(x_cls_i,size(x_2d))); title('modified Newtons');
    case 4
        [x_cls_i] = opt_bfgs(obj,x0,opt,y,h);
        figure, imshow(reshape(x_cls_i,size(x_2d))); title('BFGS');
    case 5
        [x_cls_i] = opt_lbfgs(obj,x0,opt,y,h);
        figure, imshow(reshape(x_cls_i,size(x_2d))); title('L-BFGS');
    case 6
        [x_cls_i] = opt_cg(obj,x0,opt,y,h);
        figure, imshow(reshape(x_cls_i,size(x_2d))); title('conjugate gradient');
end
psnr_cls_i = psnr(x_2d,x_cls_i,1);
fprintf('psnr: %.2f\n',psnr_cls_i);
