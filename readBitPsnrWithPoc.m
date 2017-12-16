function [ matBitPsnr ] = readBitPsnrWithPoc( dirLogFile, listPoc, isLibPic )
%Author: ylonge.
%Function: read bit and PSNR of picture in given poc range from log file.
%   --dirLogFile: directory of log file.
%   --listPoc: poc list of picture.
%   --isLibPic: true for library pictures and false for sequence pictures.
%   --matBitPsnr: N*4 matrix to store the data read, each column for dataPoc dataType dataQP dataBit dataPSNRY dataPSNRU dataPSNRV
fidLogFile = fopen(dirLogFile, 'r');
if(fidLogFile == -1)
    fprintf(2, ferror(fidLogFile));
    return;
end
numPoc = length(listPoc);
beginSymbolPoc = 'C';
endSymbolPoc = 'T';
beginSymbolType = '('; % only read one letter.
beginSymbolQP = 'nQP'; offsetQP = 9;% and offset 9 to QP.
endSymbolQP = ')';
beginSymbolBit = ')';
endSymbolBit = 'b';
beginSymbolPSNRY = 'Y';
beginSymbolPSNRU = 'U';
beginSymbolPSNRV = 'V';
beginSymbolLibLog = '####### log of library picture #######';
beginSymbolSeqLog = '####### log of sequence picture #######';
numSymbolPSNR = 8;
countDataWanted = 0;
matBitPsnr = [];
skipFlag = false;
while(~feof(fidLogFile))
    strLineExtract = fgetl(fidLogFile);

    if(isempty(strLineExtract))
        continue;
    end
    % case 1: early terminate when reading library log and finding sequence log.
    if isLibPic && strcmp(strLineExtract, beginSymbolSeqLog)
        break;
    end
    % case 2: skip when reading sequence log and finding library log.
    if ~isLibPic && strcmp(strLineExtract, beginSymbolLibLog)
        skipFlag = true;
        continue;
    end
    % case 3: do not skip when reading sequence log and finding sequence log.
    if skipFlag && ~isLibPic && strcmp(strLineExtract, beginSymbolSeqLog)
        skipFlag = false;
        continue;
    end

    if skipFlag
        continue;
    end
    
    idxBeginPoc = find(strLineExtract == beginSymbolPoc, 1, 'first') + 1;
    idxEndPoc = find(strLineExtract == endSymbolPoc, 1, 'first') - 1;
    idxBeginType = find(strLineExtract == beginSymbolType, 1) + 2;
    idxBeginQP = strfind(strLineExtract, beginSymbolQP) + offsetQP;
    idxEndQP = find(strLineExtract == endSymbolQP, 1) - 1;
    idxBeginBit = find(strLineExtract == beginSymbolBit) + 1;
    idxEndBit = find(strLineExtract == endSymbolBit) - 1;
    idxBeginPSNRY = find(strLineExtract == beginSymbolPSNRY) + 1;
    idxBeginPSNRU = find(strLineExtract == beginSymbolPSNRU) + 1;
    idxBeginPSNRV = find(strLineExtract == beginSymbolPSNRV) + 1;
    if(isempty(idxBeginPoc) || isempty(idxBeginBit) || isempty(idxEndBit) || isempty(idxBeginPSNRY) || isempty(idxBeginPSNRU) || isempty(idxBeginPSNRV))
        continue;
    end
    dataPoc = str2num(strLineExtract(idxBeginPoc: idxEndPoc));
    % check whether the poc is in the list.
    if isempty(find(listPoc == dataPoc, 1))
        continue;
    end
    
    strType = strLineExtract(idxBeginType);
    switch strType
        case 'I'
            dataType = 0;
        case 'P'
            dataType = 1;
        case 'B'
            dataType = 2;
        otherwise
            error('reading wrong slice type.');
    end
    dataQP = str2num(strLineExtract(idxBeginQP: idxEndQP));
    dataBit = str2num(strLineExtract(idxBeginBit: idxEndBit));
    dataPSNRY = str2num(strLineExtract(idxBeginPSNRY: idxBeginPSNRY + numSymbolPSNR - 1));
    dataPSNRU = str2num(strLineExtract(idxBeginPSNRU: idxBeginPSNRU + numSymbolPSNR - 1));
    dataPSNRV = str2num(strLineExtract(idxBeginPSNRV: idxBeginPSNRV + numSymbolPSNR - 1));
    matBitPsnr = [matBitPsnr; dataPoc dataType dataQP dataBit dataPSNRY dataPSNRU dataPSNRV];
    countDataWanted = countDataWanted + 1;
    if(countDataWanted == numPoc && numPoc ~= 0)
        break;
    end
end
fclose(fidLogFile);
end