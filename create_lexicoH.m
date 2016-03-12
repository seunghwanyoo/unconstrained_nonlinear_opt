function [H,h_2d] = create_lexicoH(x_2d,h0_2d)
% Given 2D blur kernel, create the linear operator matrix H
% Input:
%   x_2d: 2D image
%   h0_2d: 2D (blur) kernel
% Output:
%   H: operation matrix H with lexicographic notation
%   h_2d: 2D kernel with the size of x_2d
%
% Author: Seunghwan Yoo

disp('create_lexicoH()');
[m,n] = size(x_2d);
h_2d = create_h2d(x_2d,h0_2d);

h_row = zeros(m,m*n);
for i=0:m-1
    h_2d_shift = circshift(h_2d,[i,0]);
    h_row(i+1,:) = h_2d_shift(:)';
end
for j=0:n-1
    h_row_shift = circshift(h_row,[0,j*m]);
    H(j*m+1:(j+1)*m,:) = h_row_shift;
end