function [ listAverBitPsnr, listAverBitPsnrKeyFrames, listTotalSSE, listTotalSSEKeyFrames, listNumFrames ] = readBitPsnrSingleLib( dirLogFile, listOrgLib, listKeyPic, listRefLib, frameRate, flagLib )
%Author: ylonge.
%Function: Extract average bit and PSNR of frames that reference library pictures, respectively.
%          Note that the leading pictures of key pictures do not only reference the key picture but also reference pictures in the former GOP.
%          In this function, leading pictures are gathered with the former key picture in the encoding order.
%   --dirLogFile: directory of log file.
%   --listOrgLib: list carring original poc of library pictures. 
%   --listKeyPic: list carring poc of key pictures.
%   --listRefLib: list carring library index for each key picture.
%   --listAverBitPsnr: N*5 matrix to store the data extracted, each column for bit, PSNR-Y, PSNR-U, PSNR-V, PSNR-YUV.
%   --listAverBitPsnrKeyFrames: N*5 matrix to store only the key frame data extracted, each column for bit, PSNR-Y, PSNR-U, PSNR-V, PSNR-YUV.
%   --listTotalSSE: the summary squared error of all frames.
%   --listTotalSSEKeyFrames: the summary squared error of key frames.

% prepare library information.
numLibPic = length(listOrgLib);
numKeyPic = length(listKeyPic);

% read bit and psnr of each picture.
[listBitPsnrAllPic, listPocAll] = readLineBitPsnr( dirLogFile );
if flagLib
    listBitPsnrLibPic = listBitPsnrAllPic(1: numLibPic, :);
    listBitPsnrSeqPic = listBitPsnrAllPic((numLibPic + 1): end, :);
    listPoc = listPocAll((numLibPic + 1): end);
else
    listBitPsnrLibPic = 0;
    listBitPsnrSeqPic = listBitPsnrAllPic;
    listPoc = listPocAll;
end

% compute yuv psnr.
maxVar = 2^8 - 1;
weight = [4;1;1]; % for 420 component type.
sumMse = zeros(numLibPic, 1);

% gather the pictures referencing one library picture together.
listSumBitPsnr = zeros(numLibPic, 4);
listCountFrames = zeros(numLibPic, 1);
listSumBitPsnrKeyFrames = zeros(numLibPic, 4);
listFreq = zeros(numLibPic, 1);
listTotalSSE = zeros(numLibPic, 3);
listTotalSSEKeyFrames = zeros(numLibPic, 3);

for idxKeyPic = 1: numKeyPic
    idxLibPic = listRefLib(idxKeyPic) + 1;
    startPos = find(listPoc == listKeyPic(idxKeyPic));
    if idxKeyPic < numKeyPic
        endPos = find(listPoc == listKeyPic(idxKeyPic + 1)) - 1;
    else
        endPos = length(listPoc);
    end
    if length(startPos) ~= 1 || length(endPos) ~= 1
        error('there are 2 key pictures with the same poc\n');
    end

    % collect data.
    listCountFrames(idxLibPic) = listCountFrames(idxLibPic) + endPos - startPos + 1;
    listSumBitPsnr(idxLibPic, :) = listSumBitPsnr(idxLibPic, :) + sum(listBitPsnrSeqPic(startPos: endPos, :), 1);
    listSumBitPsnrKeyFrames(idxLibPic, :) = listSumBitPsnrKeyFrames(idxLibPic, :) + listBitPsnrSeqPic(startPos, :);
    listFreq(idxLibPic) = listFreq(idxLibPic) + 1;

    % compute yuv psnr.
    mse = maxVar * maxVar ./ (10 .^ (listBitPsnrSeqPic(startPos: endPos, 2:end) ./ 10));
    averMse = mse * weight / sum(weight);
    sumMse(idxLibPic) = sumMse(idxLibPic) + sum(averMse);

    % collect sse data.
    listTotalSSE(idxLibPic, :) = listTotalSSE(idxLibPic, :) + sum(mse, 1);
    listTotalSSEKeyFrames(idxLibPic, :) = listTotalSSEKeyFrames(idxLibPic, :) + mse(1, :);
end
% add the bits of library pictures.
listSumBitPsnr(:, 1) = (listSumBitPsnr(:, 1) + listBitPsnrLibPic(:, 1)) * frameRate / 1000;
listSumBitPsnrKeyFrames(:, 1) = (listSumBitPsnrKeyFrames(:, 1) + listBitPsnrLibPic(:, 1)) * frameRate / 1000;
% compute averange bits.
listAverBitPsnr = listSumBitPsnr ./ (listCountFrames * ones(1, size(listSumBitPsnr, 2)));
listAverBitPsnrKeyFrames = listSumBitPsnrKeyFrames ./ (listFreq * ones(1, size(listSumBitPsnrKeyFrames, 2)));
% compute yuv PSNR.
yuvPsnr = 10 * log10(maxVar * maxVar ./ (sumMse ./ listCountFrames));
listAverBitPsnr = [listAverBitPsnr yuvPsnr];
listNumFrames = listCountFrames;
end

function [ matrixData, listPoc ] = readLineBitPsnr( dirLogFile )
%Author: ylonge.
%Function: Extract bit and PSNR of each frame in log file.
%   --dirLogFile: directory of log file.
%   --numDataWanted: given number of data wanted.
%   --matrixData: N*4 matrix to store the data extracted, each column for bit, PSNR-Y, PSNR-U, PSNR-V. NOTE: there is no yuv psnr for each picture.
fidLogFile = fopen(dirLogFile, 'r');
if(fidLogFile == -1)
    fprintf(2, ferror(fidLogFile));
    return;
end
symbol = 'POC';
matrixData = [];
listPoc = [];
while(~feof(fidLogFile))
    strLineExtract = fgetl(fidLogFile);
    if(~isempty(strfind(strLineExtract, symbol)))
        lineVal = sscanf(strLineExtract, 'POC %f %*[^)]) %f bits [Y %f dB U %f dB V %f dB]');
        bitPsnr = lineVal(2: end)'; % the first two values are number and string.
        matrixData = [matrixData; bitPsnr];
        listPoc = [listPoc; lineVal(1)];
    end
end
fclose(fidLogFile);
end