function [ listSeq ] = readSeqList( dirFile )
%Author: ylonge.
%Date: 2017/11/16
%Function: Read sequence list from file.
%   --dirFile: directory of file.
%   --listSeq: cell containing name of sequences.

fidFile = fopen(dirFile, 'r');
if(fidFile == -1)
    fprintf(2, ferror(fidFile));
    return;
end

listSeq = {};
while(~feof(fidFile))
    strLineExtract = fgetl(fidFile);
    if(isempty(strLineExtract) || ~isempty(strfind(strLineExtract, '#')))
        continue;
    end

    %% text type: seqName intraPeriod frameNum seqPath (can be updated.)
    seq = sscanf(strLineExtract, '%s%*[\n]');
    num = sscanf(strLineExtract, '%*s%d%d');
    path = sscanf(strLineExtract, '%*s%*s%*s%s');

    % read resolution from sequence name.
    symbol = '_';
    symbolRes = 'x';
    midPos = find(seq == symbolRes);
    if length(midPos) > 1
        err('yuv file name has more than one x and fail to read resolution.\n');
    end
    startPos = find(seq(1: midPos) == symbol, 1, 'last');
    endPos = find(seq(midPos: end) == symbol, 1, 'first') + midPos - 1;
    width = str2num(seq((startPos + 1): (midPos - 1)));
    height = str2num(seq((midPos + 1): (endPos - 1)));
    
    % exception for 720p sequences.
    if(strfind(seq,'720p'))
        width = 1280;
        height = 720;
    end

    listSeq = [listSeq; {path} {seq} {num(1)} {num(2)} {width} {height}];
end
fclose(fidFile);
end