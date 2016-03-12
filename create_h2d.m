function [h_2d] = create_h2d(x_2d,h0_2d)
% Given small 2D blur kernel, expand to the image size
% Input:
%   x_2d: 2D image
%   h0_2d: 2D (blur) kernel
% Output:
%   h_2d: 2D kernel with the size of x_2d
%
% Author: Seunghwan Yoo

[m,n] = size(x_2d);
h1 = zeros(m,n);

[mh,nh] = size(h0_2d);
h1(1:mh,1:nh) = h0_2d; % zero-padded blur kernel
h_2d = circshift(h1,[-floor(mh/2),-floor(nh/2)]);