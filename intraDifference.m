function [ errAver, errTotal, errMat ] = intraDifference( trgPic, blockSize )
%Author: ylonge.
%Function: This function is used to process intra difference.
%   --trgPic: sample of target picture. only luma is used in common situation.
%   --blockSize: size of block based on which the intra difference is processed.
%   --err: the difference between two pictures after motion compensation.

if blockSize <= 0
    return;
end

[width, height] = size(trgPic);
numBlockInWidth = ceil(width / blockSize);
numBlockInHeight = ceil(height / blockSize);

errTotal = 0;
errMat = zeros(numBlockInWidth, numBlockInHeight);

for idxBlockInWidth = 1: numBlockInWidth
    startTrgPosW = blockSize * (idxBlockInWidth - 1) + 1;
    endTrgPosW = blockSize * idxBlockInWidth;
    intervalTrgWidth = startTrgPosW: endTrgPosW;

    for idxBlockInHeight = 1: numBlockInHeight
        startTrgPosH = blockSize * (idxBlockInHeight - 1) + 1;
        endTrgPosH = blockSize * idxBlockInHeight;
        intervalTrgHeight = startTrgPosH: endTrgPosH;
        trgBlock = fetchBlock(trgPic, intervalTrgWidth, intervalTrgHeight);
        
        %% compute the minimal difference.
        % errMin = intraDifferencePerBlock( trgPic, trgBlock, idxBlockInWidth, idxBlockInHeight );
        errMin = intraDiffGradientPerBlock( trgPic, idxBlockInWidth, idxBlockInHeight, blockSize );

        % collect data.
        errMat(idxBlockInWidth, idxBlockInHeight) = errMin;
        errTotal = errTotal + errMin;
    end
end

errAver = errTotal / numBlockInWidth / numBlockInHeight / blockSize / blockSize;

end