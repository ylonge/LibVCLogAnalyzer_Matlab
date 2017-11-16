function [yuvPsnr] = computeYuvPsnr(yPsnr, uPsnr, vPsnr, componentType)
% Author: Ylonge.
% Date: 2017/11/3.
% Function: compute yuv PSNR based on PSNR of all compenents.
%   --componentType: string indicates the type of component, i.e. '420'.

switch componentType
	case '420'
		weight = [4; 1; 1];
	otherwise
		error('component type is not given\n');
end
numComponent = 3;
numBit = 8;
maxVal = 2 ^ numBit - 1;
listPsnr = [yPsnr; uPsnr; vPsnr];
listMse = zeros(3, 1);
for idxComp = 1: numComponent
	listMse(idxComp) = maxVal * maxVal / (10 ^ (listPsnr(idxComp) / 10));
end
yuvMse = sum(listMse .* weight);
yuvPsnr = 10 * log10(maxVal * maxVal / yuvMse);
end