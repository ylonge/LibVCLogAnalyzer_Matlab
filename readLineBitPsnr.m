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