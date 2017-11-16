function [ gainPerformance ] = bdRateComparation( matrixBdRateAnchor, matrixBdRateTest )
%Author: ylonge.
%Function: This function is used to compare performance of the
%bit-distortion rate of anchor and test.
%   --matrixBdRateAnchor: 4*5 matrix of anchor, each row is bit rate, PSNR Y, PSNR
%   U, PSNR V and PSNR YUV with different Qp ranging from 22 to 37.
%   --matrixBdRateTest: 4*5 matrix of test, each row is bit rate, PSNR Y, PSNR
%   U, PSNR V and PSNR YUV with different Qp ranging from 22 to 37.
%   --gainPerformance: comparative result of performance, using %, positive for
%   loss and negtive for gain.

rateAnchor = matrixBdRateAnchor(:, 1);
rateTest = matrixBdRateTest(:, 1);
numComponent = size(matrixBdRateAnchor, 2) - 1;
gainPerformance = zeros(1, numComponent);
for idxComponent = 1: numComponent
    distortionAnchor = matrixBdRateAnchor(:, idxComponent + 1);
    distortionTest = matrixBdRateTest(:, idxComponent + 1);
    gain = bdRateSingle(rateAnchor, distortionAnchor, rateTest, distortionTest);
    gainPerformance(idxComponent) = gain;
end

end

function [ gain ] = bdRateSingle( rateAnchor, distortionAnchor, rateTest, distortionTest )
% This child funtion is used to compute each PSNR gain under given rate.
minPSNR = max(min(distortionAnchor), min(distortionTest));
maxPSNR = min(max(distortionAnchor), max(distortionTest));

vA = bdRint(rateAnchor, distortionAnchor, minPSNR, maxPSNR);
vB = bdRint(rateTest, distortionTest, minPSNR, maxPSNR);

average = (vB - vA) / (maxPSNR - minPSNR);
gain = 10.^average - 1;
end

function [ result ] = bdRint(rate, distortion, low, high)
% This child funciton is used by function bdRateSingle();
for i = 1: 4
    logRate(i) = log10(rate(5 - i));
    logDistortion(i) = distortion(5 - i);
end

for i = 1: 3
    H(i) = logDistortion(i + 1) - logDistortion(i);
    delta(i) = (logRate(i + 1) - logRate(i)) / H(i);
end

d(1) = pchipend(H(1), H(2), delta(1), delta(2));
for i = 2: 3
    d(i) = (3 * H(i - 1) + 3 * H(i)) / ((2 * H(i) + H(i - 1)) / delta(i - 1) + (H(i) + 2 * H(i - 1)) / delta(i));
end
d(4) = pchipend(H(3), H(2), delta(3), delta(2));
for i = 1: 3
    c(i) = (3 * delta(i) - 2 * d(i) - d(i + 1)) / H(i);
    b(i) = (d(i) - 2 * delta(i) + d(i + 1)) / (H(i) * H(i));
end

result = 0;
for i = 1: 3
    s0 = logDistortion(i);
    s1 = logDistortion(i + 1);
    
    s0 = max(s0, low);
    s0 = min(s0, high);
    
    s1 = max(s1, low);
    s1 = min(s1, high);
    
    s0 = s0 - logDistortion(i);
    s1 = s1 - logDistortion(i);
    
    if(s1 > s0)
        result = result + (s1 - s0) * logRate(i);
        result = result + (s1 * s1 - s0 * s0) * d(i) / 2;
        result = result + (s1 * s1 * s1 - s0 * s0 * s0) * c(i) / 3;
        result = result + (s1 * s1 * s1 * s1 - s0 * s0 * s0 * s0) * b(i) / 4;
    end
end
end

function [ d ] = pchipend( h1, h2, delta1, delta2 )
d = ((2 * h1 + h2) * delta1 - h1 * delta2) / (h1 + h2);
if(d * delta1 < 0)
    d = 0;
else if((delta1 * delta2 < 0) && (abs(d) > abs(3 * delta1)))
        d = 3 * delta1;
    end
end
end

