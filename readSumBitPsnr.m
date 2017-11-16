function [ bitPsnr ] = readSumBitPsnr( dirLogFile )
%Author: ylonge.
%Function: Read summary of bit and PSNR in log file.
%   --dirLogFile: directory of log file.
%   --bitPsnr: 1*5 matrix to store the data, for bit, PSNR-Y, PSNR-U, PSNR-V, PSNR-YUV.

fidLogFile = fopen(dirLogFile, 'r');
if(fidLogFile == -1)
    fprintf(2, ferror(fidLogFile));
    return;
end
symbol = 'SUMMARY';

while(~feof(fidLogFile))
    strLineExtract = fgetl(fidLogFile);
    if(~isempty(strfind(strLineExtract, symbol)))
        strLineExtract = fgetl(fidLogFile);
        strLineExtract = fgetl(fidLogFile);
        lineVal = sscanf(strLineExtract, '%f%s%f%f%f%f%f');
        bitPsnr = lineVal(3: end); % the first two values are number and string.
        break;
    end
end
fclose(fidLogFile);
end