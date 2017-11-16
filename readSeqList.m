function [ listSeq ] = readSeqList( dirFile )
%Author: ylonge.
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
    seq = sscanf(strLineExtract, '%s%*[\n]');
    num = sscanf(strLineExtract, '%*s%d%d');
    path = sscanf(strLineExtract, '%*s%*s%*s%s');
    listSeq = [listSeq; {path} {seq} {num(1)} {num(2)}];
end
fclose(fidFile);
end