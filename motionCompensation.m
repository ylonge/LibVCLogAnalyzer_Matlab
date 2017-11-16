function [ errInterAver, errInterTotal, errIntraAver, errIntraTotal, percentIntra, mvCell, errMat ] = motionCompensation( refPic, trgPic, blockSize, searchRange )
%Author: ylonge.
%Function: This function is used to process motion compensation and compute the error.
%   --refPic: sample of reference picture. only luma is used in common situation.
%   --trgPic: sample of target picture. only luma is used in common situation.
%   --blockSize: size of block based on which the motion compentation is processed.
%   --searchRange: range of search.
%   --errInterTotal: the difference between two pictures after motion compensation.
%   --errIntraTotal: the average difference of intra blocks.
%	--mvCell: cell carrying motion vectors of all blocks.
%	--errMat: matrix carrying error of all blocks.

if log2(searchRange) ~= fix(log2(searchRange))
    error('search range is not exponent of 2\n');
end

[width, height] = size(trgPic);
numBlockInWidth = ceil(width / blockSize);
numBlockInHeight = ceil(height / blockSize);

errInterTotal = 0;
errIntraTotal = 0;
mvCell = cell(numBlockInWidth, numBlockInHeight);
errMat = zeros(numBlockInWidth, numBlockInHeight);
countIntraBlock = 0;

for idxBlockInWidth = 1: numBlockInWidth
    startTrgPosW = blockSize * (idxBlockInWidth - 1) + 1;
    endTrgPosW = blockSize * idxBlockInWidth;
    intervalTrgWidth = startTrgPosW: endTrgPosW;

    for idxBlockInHeight = 1: numBlockInHeight
        startTrgPosH = blockSize * (idxBlockInHeight - 1) + 1;
        endTrgPosH = blockSize * idxBlockInHeight;
        intervalTrgHeight = startTrgPosH: endTrgPosH;
        trgBlock = fetchBlock(trgPic, intervalTrgWidth, intervalTrgHeight);
        
        % search and compute error.
        SearchCenter = [startTrgPosW; startTrgPosH];
        [errInter, mv] = motionCompensationBlock(refPic, trgBlock, searchRange, SearchCenter);

        % decide which mode is used.
        % errIntra = intraDifferencePerBlock( trgPic, trgBlock, idxBlockInWidth, idxBlockInHeight );
        errIntra = intraDiffGradientPerBlock( trgPic, idxBlockInWidth, idxBlockInHeight, blockSize );
        if errIntra < errInter
        	countIntraBlock = countIntraBlock + 1;
        	errIntraTotal = errIntraTotal + errIntra;
        	errMat(idxBlockInWidth, idxBlockInHeight) = errIntra;
        else
        	errInterTotal = errInterTotal + errInter;
        	mvCell(idxBlockInWidth, idxBlockInHeight) = {mv};
        	errMat(idxBlockInWidth, idxBlockInHeight) = errInter;
        end
    end
end

errInterAver = errInterTotal / (numBlockInWidth * numBlockInHeight - countIntraBlock) / blockSize / blockSize;
if countIntraBlock >= 0
	errIntraAver = errIntraTotal / countIntraBlock / blockSize / blockSize;
end

percentIntra = countIntraBlock / (numBlockInWidth * numBlockInHeight);

end

function [err, mv] = motionCompensationBlock(refPic, trgBlock, searchRange, initSearchCenter)

numIter = log2(searchRange) + 2; % for position of 0, 2^0, 2^1, 2^2......

SearchCenterW = initSearchCenter(1);
SearchCenterH = initSearchCenter(2);
blockSize = size(trgBlock, 1);

stepThr = 3;
rangeAll = 3;

mvBest = zeros(2, 1);
errMin = Inf;

while (1)
    mvCur = zeros(2, 1);
    idxIterCur = 0;
    for idxIter = 1: numIter
        numIterW = 1;
        numIterH = 1;
        listPosW = 0;
        listPosH = 0;
        if idxIter > 1
            numIterW = 3;
            numIterH = 3;
            listPosW = [-2^(idxIter - 2); 0; 2^(idxIter - 2)];
            listPosH = [-2^(idxIter - 2); 0; 2^(idxIter - 2)];
        end

        for idxIterW = 1: numIterW
            posW = listPosW(idxIterW) + SearchCenterW;
            intervalW = posW: (posW + blockSize - 1);
            for idxIterH = 1: numIterH
                posH = listPosH(idxIterH) + SearchCenterH;
                intervalH = posH: (posH + blockSize - 1);

                refBlock = fetchBlock(refPic, intervalW, intervalH);
                errCur = diffBlock(refBlock, trgBlock);
                if errCur < errMin
                    errMin = errCur;
                    mvCur(1) = listPosW(idxIterW);
                    mvCur(2) = listPosH(idxIterH);
                    idxIterCur = idxIter;
                end
            end
        end
    end

    % check if the iteration step is less than 1.
    if idxIterCur == 2
        if abs(mvCur(1)) > 2 || abs(mvCur(2)) > 2
            error('ERROR: the search process has something wrong.\n');
        end
        mvCheck1 = mvCur;
        mvCheck2 = mvCur;
        listPos = [1; 2];

        zeroPos = listPos(mvCur == 0);
        if isempty(zeroPos)
            negPos = listPos(mvCur < 0);
            if isempty(negPos)
                mvCheck1(1) = mvCheck1(1) + 1;
                mvCheck2(2) = mvCheck2(2) + 1;
            elseif length(negPos) == 1
                posPos = 3 - negPos;
                mvCheck1(negPos) = mvCheck1(negPos) - 1;
                mvCheck2(posPos) = mvCheck2(posPos) + 1;
            else
                mvCheck1(1) = mvCheck1(1) - 1;
                mvCheck2(2) = mvCheck2(2) - 1;
            end
        elseif length(zeroPos) == 1
            mvCheck1(zeroPos) = 1;
            mvCheck2(zeroPos) = -1;
            nonZeroPos = 3 - zeroPos;
            if(mvCheck1(nonZeroPos) >= 0)
                mvCheck1(nonZeroPos) = mvCheck1(nonZeroPos) + 1;
                mvCheck2(nonZeroPos) = mvCheck2(nonZeroPos) + 1;
            else
                mvCheck1(nonZeroPos) = mvCheck1(nonZeroPos) - 1;
                mvCheck2(nonZeroPos) = mvCheck2(nonZeroPos) - 1;
            end                
        else
            error('ERROR: the iteration 2 find the center.\n')
        end

        for idx = 1: 2
            if idx == 1
                mvCheck = mvCheck1;
            else
                mvCheck = mvCheck2;
            end

            posW = mvCheck(1) + SearchCenterW;
            posH = mvCheck(2) + SearchCenterH;
            intervalW = posW: (posW + blockSize - 1);
            intervalH = posH: (posH + blockSize - 1);

            refBlock = fetchBlock(refPic, intervalW, intervalH);
            errCur = diffBlock(refBlock, trgBlock);
            if errCur < errMin
                errMin = errCur;
                mvCur(1) = mvCheck(1);
                mvCur(2) = mvCheck(2);
            end
        end
    end

    % check whether the iteration is larger than threshold. full search is processed.
    if idxIterCur > stepThr
        mvCheckInit = mvCur;
        for idxW = -rangeAll: rangeAll
            for idxH = -rangeAll: rangeAll
                mvCheck = mvCheckInit + [idxW; idxH];

                posW = mvCheck(1) + SearchCenterW;
                posH = mvCheck(2) + SearchCenterH;
                intervalW = posW: (posW + blockSize - 1);
                intervalH = posH: (posH + blockSize - 1);

                refBlock = fetchBlock(refPic, intervalW, intervalH);
                errCur = diffBlock(refBlock, trgBlock);
                if errCur < errMin
                    errMin = errCur;
                    mvCur(1) = mvCheck(1);
                    mvCur(2) = mvCheck(2);
                end
            end
        end
    end

    if mvCur(1) == 0 && mvCur(2) == 0
        break;
    end

    mvBest = mvCur + mvBest;
    SearchCenterW = initSearchCenter(1) + mvBest(1);
    SearchCenterH = initSearchCenter(2) + mvBest(2);
end
err = errMin;
mv = mvBest;
end