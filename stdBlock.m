function [stdBlk] = stdBlock(pic, sizeBlock)
[wp, hp] = size(pic);
wb = wp / sizeBlock;
hb = hp / sizeBlock;
sumStdBlk = 0;
for iw = 1: wb
    sw = (iw - 1) * sizeBlock + 1;
    ew = iw * sizeBlock;
    rw = sw: ew;
    for ih = 1: hb
        sh = (ih - 1) * sizeBlock + 1;
        eh = ih * sizeBlock;
        rh = sh: eh;

        block = fetchBlock(pic, rw, rh);
        stdEachBlock = sqrt(mean((block(:) - mean(block(:))) .^ 2));
        sumStdBlk = sumStdBlk + stdEachBlock;
    end
end
stdBlk = sumStdBlk / wb / hb;
end

function [refBlock] = fetchBlock(refPic, intervalW, intervalH)

wb = length(intervalW);
hb = length(intervalH);
refBlock = zeros(wb, hb) + 128;
for idx = 1: wb
    if intervalW(idx) <= 0
        continue;
    end
    maskFill = intervalH > 0;
    sampleFill = refPic(intervalH(maskFill));
    rowFill = intervalW(idx);
    refBlock(rowFill, maskFill) = sampleFill;
end
end