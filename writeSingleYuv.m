function [ y,u,v ] = writeSingleYuv( dirFile, y, u, v, bits )
%Author: ylonge.
%Function: Read picture yuv from yuv file.
%   --dirLogFile: directory of yuv file.
%   --poc: the poc of wanted picture.

fidFile = fopen(dirFile, 'wb');
if(fidFile == -1)
    fprintf(2, ferror(fidFile));
    return;
end

width = size(y, 1);
height = size(y, 2);

if bits <= 8
    lengthY = height * width;
    line = reshape(y, lengthY, 1);
    fwrite(fidFile, line);

    lengthUV = height * width / 4;
    line = reshape(u, lengthUV, 1);
    fwrite(fidFile, line);

    line = reshape(v, lengthUV, 1);
    fwrite(fidFile, line);
elseif bits > 8
    lengthY = height * width * 2;
    lineComb = reshape(y, width * height, 1);
    lineLow = rem(lineComb, 256);
    lineHigh = fix(lineComb./256);
    line = zeros(lengthY, 1);
    line(1:2:end) = lineLow;
    line(2:2:end) = lineHigh;
    fwrite(fidFile, line);

    lengthUV = lengthY / 4;
    lineComb = reshape(u, width * height / 4, 1);
    lineLow = rem(lineComb, 256);
    lineHigh = fix(lineComb./256);
    line = zeros(lengthUV, 1);
    line(1:2:end) = lineLow;
    line(2:2:end) = lineHigh;
    fwrite(fidFile, line);

    lineComb = reshape(v, width * height / 4, 1);
    lineLow = rem(lineComb, 256);
    lineHigh = fix(lineComb./256);
    line = zeros(lengthUV, 1);
    line(1:2:end) = lineLow;
    line(2:2:end) = lineHigh;
    fwrite(fidFile, line);
end

fclose(fidFile);
end