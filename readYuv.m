function [ y,u,v ] = readYuv( dirFile, poc, width, height )
%Author: ylonge.
%Function: Read picture yuv from yuv file.
%   --dirLogFile: directory of yuv file.
%   --poc: the poc of wanted picture.

fidFile = fopen(dirFile, 'r');
if(fidFile == -1)
    fprintf(2, ferror(fidFile));
    return;
end

if nargin <= 2
    % read resolution from dirFile.
    symbol = '_';
    symbolRes = 'x';
    midPos = find(dirFile == symbolRes);
    if length(midPos) > 1
        err('yuv file name has more than one x and fail to read resolution.\n');
    end
    startPos = find(dirFile(1: midPos) == symbol, 1, 'last');
    endPos = find(dirFile(midPos: end) == symbol, 1, 'first') + midPos - 1;
    width = str2num(dirFile((startPos + 1): (midPos - 1)));
    height = str2num(dirFile((midPos + 1): (endPos - 1)));
    
    % exception for 720p sequences.
    if(strfind(dirFile,'720p'))
        width = 1280;
        height = 720;
    end
end
if nargin <= 1
    err('poc is not given!\n');
end

skip = poc * height * width * 1.5;
fseek(fidFile, skip, 'bof');

lengthY = height * width;
line = fread(fidFile, lengthY);
y = reshape(line, width, height);

lengthUV = height * width / 4;
line = fread(fidFile, lengthUV);
u = reshape(line, width/2, height/2);

line = fread(fidFile, lengthUV);
v = reshape(line, width/2, height/2);


fclose(fidFile);
end