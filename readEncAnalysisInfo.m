%% read library information file.
function [ encAnalInfo ] = readEncAnalysisInfo( encAnalInfoFile )
% Author: Ylonge.
% Date: 2018/1/16.
% Function: Read encode analysis information file generated by LibVC encoder based on HM-16-15.
%   --encAnalInfo: struct storing analysis information.
%   --encAnalInfoFile: file containing analysis information.

fidEncAnalInfoFile = fopen( encAnalInfoFile );
if(fidEncAnalInfoFile == -1)
    fprintf(2, ferror(fidEncAnalInfoFile));
    return;
end
encAnalInfo = cell(2, 1); % 1 for library pictures and 2 for key pictures.
keyAnalInfo = [];
libAnalInfo = [];

while(~feof(fidEncAnalInfoFile))
    strLineExtract = fgetl(fidEncAnalInfoFile);
    if(isempty(strLineExtract))
        continue;
    end
    
    if (~isempty(strfind(strLineExtract, 'keyG')))
        lineVal = sscanf(strLineExtract(5:end), '%f%f%f%f%f%f%f%f%f');
        keyAnalInfo = [keyAnalInfo; lineVal'];
    elseif (~isempty(strfind(strLineExtract, 'libG')))
        lineVal = sscanf(strLineExtract(5:end), '%f%f%f%f');
        libAnalInfo = [libAnalInfo; lineVal'];
    end
end
encAnalInfo{1} = libAnalInfo;
encAnalInfo{2} = keyAnalInfo;

fclose(fidEncAnalInfoFile);
end