function [ listAverBitPsnr, listAverBitPsnrKeyFrames, scalePicDivKeyRD, listNumFrames, listBitDistAllPicAllLib ] = readBitPsnrSingleLib( dirLogFile, listOrgLib, listKeyPic, listRefLib, frameRate, flagLib )
%Author: ylonge.
%Function: Extract average bit and PSNR of frames that reference library pictures, respectively.
%          Note that the leading pictures of key pictures do not only reference the key picture but also reference pictures in the former GOP.
%          In this function, leading pictures are gathered with the former key picture in the encoding order.
%Input:
%   --dirLogFile: directory of log file.
%   --listOrgLib: list carring original poc of library pictures. 
%   --listKeyPic: list carring poc of key pictures.
%   --listRefLib: list carring library index for each key picture.
%Output:
%   --listAverBitPsnr: N*5 matrix to store the data extracted, each column for bit, PSNR-Y, PSNR-U, PSNR-V, PSNR-YUV.
%   --listAverBitPsnrKeyFrames: N*5 matrix to store only the key frame data extracted, each column for bit, PSNR-Y, PSNR-U, PSNR-V, PSNR-YUV.
%   --scalePicDivKeyRD: N*2 matrix carrying scale of rate and distortion of normal pictures in one intra period divide those of key pictures.
%   --listNumFrames: list of number of frames for each library picture.

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

% gather the pictures referencing one library picture together.
listSumBitPsnr = zeros(numLibPic, 4);
listCountFrames = zeros(numLibPic, 1);
listSumBitPsnrKeyFrames = zeros(numLibPic, 4);
listFreq = zeros(numLibPic, 1);
scalePicDivKeyRD = zeros(numKeyPic, 4);
listBitDistAllPicAllLib = cell(numLibPic, 1);
for i = 1: numLibPic
    numKeyPicPerLib = sum(listRefLib == (i - 1));
    listBitDistAllPicAllLib{i} = cell(numKeyPicPerLib, 1);
end

countPerLib = zeros(numLibPic, 1);
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
    countPerLib(idxLibPic) = countPerLib(idxLibPic) + 1;

    % collect data.
    listCountFrames(idxLibPic) = listCountFrames(idxLibPic) + endPos - startPos + 1;
    listSumBitPsnr(idxLibPic, :) = listSumBitPsnr(idxLibPic, :) + sum(listBitPsnrSeqPic(startPos: endPos, 1:4), 1);
    listSumBitPsnrKeyFrames(idxLibPic, :) = listSumBitPsnrKeyFrames(idxLibPic, :) + listBitPsnrSeqPic(startPos, 1:4);
    listFreq(idxLibPic) = listFreq(idxLibPic) + 1;
    listBitDistAllPicAllLib{idxLibPic}{countPerLib(idxLibPic)} = listBitPsnrSeqPic(startPos: endPos, :);

    % collect scale.
    tempList = [listBitPsnrSeqPic(startPos: endPos, 1) psnr2mse(listBitPsnrSeqPic(startPos: endPos, 2:4))];
    scalePicDivKeyRD(idxKeyPic, :) = sum(tempList, 1) ./ (endPos - startPos + 1) ./ tempList(1, :);
end
% add the bits of library pictures.
listSumBitPsnr(:, 1) = (listSumBitPsnr(:, 1) + listBitPsnrLibPic(:, 1)) * frameRate / 1000;
listSumBitPsnrKeyFrames(:, 1) = (listSumBitPsnrKeyFrames(:, 1) + listBitPsnrLibPic(:, 1)) * frameRate / 1000;
% compute averange bits.
listAverBitPsnr = listSumBitPsnr ./ (listCountFrames * ones(1, size(listSumBitPsnr, 2)));
listAverBitPsnrKeyFrames = listSumBitPsnrKeyFrames ./ (listFreq * ones(1, size(listSumBitPsnrKeyFrames, 2)));
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
        lineVal = sscanf(strLineExtract, 'POC %f %*[^Q]QP %f %*[^)]) %f bits [Y %f dB U %f dB V %f dB]');
        bitPsnr = lineVal(3: end)'; % the first two values are number and string.
        matrixData = [matrixData; bitPsnr lineVal(2)];
        listPoc = [listPoc; lineVal(1)];
    end
end
fclose(fidLogFile);
end