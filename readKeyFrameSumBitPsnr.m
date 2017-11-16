function [ bitPsnr ] = readKeyFrameSumBitPsnr( dirLogFile )
%Author: ylonge.
%Function: Read summary of bit and PSNR in log file.
%   --dirLogFile: directory of log file.
%   --bitPsnr: 1*5 matrix to store the data, for bit, PSNR-Y, PSNR-U, PSNR-V, PSNR-YUV.

fidLogFile = fopen(dirLogFile, 'r');
if(fidLogFile == -1)
    fprintf(2, ferror(fidLogFile));
    return;
end

symbolIntraPic = 'I Slices';
numIntraPic = 0;
bitPsnrIntraFrame = zeros(1, 5);

% this part is used for key frames in LibVC
symbolKeyFrame = 'P Slices';
numKeyFrame = 0;
bitPsnrKeyFrame = zeros(1, 5);

symbolLibFrame = 'Lib Pictures';
numLibFrame = 0;
bitLib =0;

while(~feof(fidLogFile))
    strLineExtract = fgetl(fidLogFile);
    if(~isempty(strfind(strLineExtract, symbolIntraPic)))
        strLineExtract = fgetl(fidLogFile);
        strLineExtract = fgetl(fidLogFile);
        lineVal = sscanf(strLineExtract, '%f%s%f%f%f%f%f');
        bitPsnrIntraFrame = lineVal(3: end); % the first two values are number and string.
        numIntraPic = lineVal(1);
    elseif(~isempty(strfind(strLineExtract, symbolKeyFrame)))
        strLineExtract = fgetl(fidLogFile);
        strLineExtract = fgetl(fidLogFile);
        lineVal = sscanf(strLineExtract, '%f%s%f%f%f%f%f');
        bitPsnrKeyFrame = lineVal(3: end); % the first two values are number and string.
        numKeyFrame = lineVal(1);
    elseif(~isempty(strfind(strLineExtract, symbolLibFrame)))
        strLineExtract = fgetl(fidLogFile);
        strLineExtract = fgetl(fidLogFile);
        lineVal = sscanf(strLineExtract, '%f%s%f%f%f%f%f');
        bitLib = lineVal(3); % the first two values are number and string.
        numLibFrame = lineVal(1);
    end
end
bitPsnr = (bitPsnrIntraFrame * numIntraPic + bitPsnrKeyFrame * numKeyFrame) / (numIntraPic + numKeyFrame);
bitPsnr(1) = bitPsnr(1) + bitLib * numLibFrame / (numIntraPic + numKeyFrame);
fclose(fidLogFile);
end