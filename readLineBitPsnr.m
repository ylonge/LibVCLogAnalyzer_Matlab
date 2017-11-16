function [ matrixData, listPoc ] = readLineBitPsnr( dirLogFile, numDataWanted )
%Author: ylonge.
%Function: Extract bit and PSNR of each frame in log file.
%   --dirLogFile: directory of log file.
%   --numDataWanted: given number of data wanted.
%   --matrixData: N*4 matrix to store the data extracted, each column for
%       bit, PSNR-Y, PSNR-U, PSNR-V.
fidLogFile = fopen(dirLogFile, 'r');
if(fidLogFile == -1)
    fprintf(2, ferror(fidLogFile));
    return;
end
beginSymbolPoc = 'C';
endSymbolPoc = 'T';
beginSymbolBit = ')';
endSymbolBit = 'b';
beginSymbolPSNRY = 'Y';
beginSymbolPSNRU = 'U';
beginSymbolPSNRV = 'V';
numSymbolPSNR = 8;
countDataWanted = 0;
matrixData = [];
listPoc = [];
while(~feof(fidLogFile))
    strLineExtract = fgetl(fidLogFile);
    if(isempty(strLineExtract))
        continue;
    end
    idxBeginPoc = find(strLineExtract == beginSymbolPoc, 1, 'first') + 1;
    idxEndPoc = find(strLineExtract == endSymbolPoc, 1, 'first') - 1;
    idxBeginBit = find(strLineExtract == beginSymbolBit) + 1;
    idxEndBit = find(strLineExtract == endSymbolBit) - 1;
    idxBeginPSNRY = find(strLineExtract == beginSymbolPSNRY) + 1;
    idxBeginPSNRU = find(strLineExtract == beginSymbolPSNRU) + 1;
    idxBeginPSNRV = find(strLineExtract == beginSymbolPSNRV) + 1;
    if(isempty(idxBeginPoc) || isempty(idxBeginBit) || isempty(idxEndBit) || isempty(idxBeginPSNRY) || isempty(idxBeginPSNRU) || isempty(idxBeginPSNRV))
        continue;
    end
    dataPoc = str2num(strLineExtract(idxBeginPoc: idxEndPoc));
    listPoc = [listPoc; dataPoc];
    
    dataBit = str2num(strLineExtract(idxBeginBit: idxEndBit));
    dataPSNRY = str2num(strLineExtract(idxBeginPSNRY: idxBeginPSNRY + numSymbolPSNR - 1));
    dataPSNRU = str2num(strLineExtract(idxBeginPSNRU: idxBeginPSNRU + numSymbolPSNR - 1));
    dataPSNRV = str2num(strLineExtract(idxBeginPSNRV: idxBeginPSNRV + numSymbolPSNR - 1));
    matrixData = [matrixData; dataBit dataPSNRY dataPSNRU dataPSNRV];
    countDataWanted = countDataWanted + 1;
    if(countDataWanted == numDataWanted && numDataWanted ~= 0)
        break;
    end
end
fclose(fidLogFile);
end