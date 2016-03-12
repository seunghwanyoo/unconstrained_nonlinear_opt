function out = psnr(img1,img2,vmax)

if nargin == 2
    if max(max(img1)) > 10
        vmax = 255;
    else
        vmax = 1;
    end
end

mse = sum((img1(:)-img2(:)).^2) / numel(img1);
out = 10 * log10(vmax*vmax/mse);