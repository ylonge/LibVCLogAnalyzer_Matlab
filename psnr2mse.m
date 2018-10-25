function [mse] = psnr2mse(psnr)
maxVar = 255;
mse = maxVar * maxVar ./ (10 .^ (psnr ./ 10));
end