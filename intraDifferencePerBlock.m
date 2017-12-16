function [ err ] = intraDifferencePerBlock( trgPic, trgBlock, idxBlockInWidth, idxBlockInHeight )
%Author: ylonge.
%Function: This function is used to process intra difference.
%   --trgPic: sample of target picture. only luma is used in common situation.
%   --trgBlock: block which is compute difference on.
%   --idxBlockInWidth: indicate the width idx of block in picture.
%   --idxBlockInHeight: indicate the height idx of block in picture.
%   --err: the difference in intra mode.

blockSize = size(trgBlock, 1);

[width, height] = size(trgPic);
numBlockInWidth = ceil(width / blockSize);
numBlockInHeight = ceil(height / blockSize);

       
%% compute the minimal difference.
% case 1: the block uses 128 as reference.
refBlockMean = zeros(blockSize, blockSize) + 128;
errMin = diffBlock(refBlockMean, trgBlock);

% case 2: the block uses left block as reference.
idxBlockInWidthRef = idxBlockInWidth - 1;
startTrgPosWRef = blockSize * (idxBlockInWidthRef - 1) + 1;
endTrgPosWRef = blockSize * idxBlockInWidthRef;
intervalTrgWidthRef = startTrgPosWRef: endTrgPosWRef;

idxBlockInHeightRef = idxBlockInHeight;
startTrgPosHRef = blockSize * (idxBlockInHeightRef - 1) + 1;
endTrgPosHRef = blockSize * idxBlockInHeightRef;
intervalTrgHeightRef = startTrgPosHRef: endTrgPosHRef;

refBlock = fetchBlock(trgPic, intervalTrgWidthRef, intervalTrgHeightRef);

errTmp = diffBlock(refBlock, trgBlock);
if errTmp < errMin
    errMin = errTmp;
end

% case 3: the block uses upper block as reference.
idxBlockInWidthRef = idxBlockInWidth;
startTrgPosWRef = blockSize * (idxBlockInWidthRef - 1) + 1;
endTrgPosWRef = blockSize * idxBlockInWidthRef;
intervalTrgWidthRef = startTrgPosWRef: endTrgPosWRef;

idxBlockInHeightRef = idxBlockInHeight - 1;
startTrgPosHRef = blockSize * (idxBlockInHeightRef - 1) + 1;
endTrgPosHRef = blockSize * idxBlockInHeightRef;
intervalTrgHeightRef = startTrgPosHRef: endTrgPosHRef;

refBlock = fetchBlock(trgPic, intervalTrgWidthRef, intervalTrgHeightRef);

errTmp = diffBlock(refBlock, trgBlock);
if errTmp < errMin
    errMin = errTmp;
end

% collect data.
err = errMin;

end