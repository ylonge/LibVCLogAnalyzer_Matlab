function [ err ] = intraDiffGradientPerBlock( trgPic, idxBlockInWidth, idxBlockInHeight, blockSize )
%Author: ylonge.
%Function: This function is used to process intra difference based on gradient of block.
%   --trgPic: sample of target picture. only luma is used in common situation.
%   --idxBlockInWidth: indicate the width idx of block in picture.
%   --idxBlockInHeight: indicate the height idx of block in picture.
%	--blockSize: size of block.
%   --err: the difference between two pictures after motion compensation.

%% prepare corrdinates of block samples.
startTrgPosW = blockSize * (idxBlockInWidth - 1) + 1;
endTrgPosW = blockSize * idxBlockInWidth;
intervalTrgWidth = startTrgPosW: endTrgPosW;

startTrgPosH = blockSize * (idxBlockInHeight - 1) + 1;
endTrgPosH = blockSize * idxBlockInHeight;
intervalTrgHeight = startTrgPosH: endTrgPosH;

%% compute gradient of each sample and collect histogram.
numQuatization = 8;
histGradDegree = zeros(numQuatization, 1);
stepDegree = 180 / numQuatization;
for idxW = 1: blockSize
	posWCurr = intervalTrgWidth(idxW);
	for idxH = 1: blockSize
		% fetch samples.
		posHCurr = intervalTrgHeight(idxH);
		sampleLeft = fetchSample(trgPic, posWCurr - 1, posHCurr);
		sampleRight = fetchSample(trgPic, posWCurr + 1, posHCurr);
		sampleUp = fetchSample(trgPic, posWCurr, posHCurr - 1);
		sampleDown = fetchSample(trgPic, posWCurr, posHCurr + 1);

		% compute degree of gradient.
		if sampleLeft == sampleRight
			gradDegree = 90;
		else
			gradDegree = atan((sampleDown - sampleUp) / (sampleRight - sampleLeft)) / pi * 180;
		end

		% compute histogram.
		idxHist = ceil((mod((gradDegree + stepDegree / 2 + 90), 180)) / stepDegree);
		histGradDegree(idxHist) = histGradDegree(idxHist) + 1;
	end
end

%% compute weight based on histogram.
maxStdHist = sqrt((blockSize * blockSize * (1 - 1 / numQuatization)) ^ 2 / numQuatization);
weightOnHist = 1 - std(histGradDegree, 1) / maxStdHist;

%% compute weighted MSE of intra mode.
trgBlock = fetchBlock(trgPic, intervalTrgWidth, intervalTrgHeight);
varBlock = var(trgBlock(:), 1) * blockSize * blockSize; % summary squared error.
err = varBlock * weightOnHist;

end

function [ sample ] = fetchSample(trgPic, posW, posH)
%% fetch sample from picture and check whether the pos exceeds edge of picture.
[w,h] = size(trgPic);

if posW <= 0 || posW > w || posH <= 0 || posH > h
	sample = 128;
else
	sample = trgPic(posW, posH);
end
end