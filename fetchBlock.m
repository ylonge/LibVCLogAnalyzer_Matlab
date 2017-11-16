function [refBlock] = fetchBlock(refPic, intervalW, intervalH)
[w,h] = size(refPic);
wb = length(intervalW);
hb = length(intervalH);
refBlock = zeros(wb, hb) + 128;
for idx = 1: wb
    if intervalW(idx) <= 0 || intervalW(idx) > w
        continue;
    end
    maskFill = intervalH > 0 & intervalH <= h;
    rowFill = intervalW(idx);
    sampleFill = refPic(rowFill, intervalH(maskFill));
    refBlock(idx, maskFill) = sampleFill;
end
end