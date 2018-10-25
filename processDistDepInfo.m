%% Function: Used for processing distortion dependency information
%%           by reading data, clustering data, show data.
%% Date: 2018/07/18
%% Author: Ylonge-ZJU.

%% Definition of constant values.
seqClass = ['A'; 'B'; 'C'; 'D'; 'E'; 'F'];
seqFilePreFix = 'SeqList';
analInfoFilePreFix = 'tempAnalysisInfo';
pathAnalInfo = '.\copy\';
seqFileSubFix = '.txt';
qpList = [19; 22; 27; 32; 37];
sizeClass = length(seqClass);
sizeBlockType = 4;
sizePicturePerSeq = 16;
sizeQP = length(qpList);

%% Pre-process to read sequence list.
seqList = cell(sizeClass, 1);
for idxClass = 1: sizeClass
    % prepare values.
	tmpFileClass = [seqFilePreFix '_' seqClass(idxClass) seqFileSubFix];
	tmpFpClass = fopen(tmpFileClass, 'r');
	% prepare memory.
    tmpSeqListPerClass = [];
	
	% read file in loop.
	while ~feof(tmpFpClass)
        tmpStrLine = fgetl(tmpFpClass);
	    tmpStrSeq = sscanf(tmpStrLine, '%s\t%*[^\n\r]');
        tmpSeqListPerClass = [tmpSeqListPerClass; {tmpStrSeq}];
	end
    
    % collect memory.
    seqList{idxClass} = tmpSeqListPerClass;

    % clear memory.
    clear tmp*;
end

%% First Step to read data with format of [depth, distEnc, distRef0, distRef1]
analInfo = cell(sizeClass, 1);
for idxClass = 1: sizeClass
	% prepare values.
	tmpSizeSeqPerClass = length(seqList{idxClass});
	% prepare memory.
	tmpAnalInfoPerClass = cell(tmpSizeSeqPerClass, sizeQP);

	% collect data in sequence loop.
	for idxSeq = 1: tmpSizeSeqPerClass

		% collect data in QP loop.
		for idxQP = 1: sizeQP
			% prepare values.
			tmpFileAnalInfoPerSeq = [pathAnalInfo seqList{idxClass}{idxSeq} '\' analInfoFilePreFix num2str(qpList(idxQP)) seqFileSubFix];
	        tmpFpAnalInfoPerSeq = fopen(tmpFileAnalInfoPerSeq, 'r');
			% prepare memory.
			tmpAnalInfoPerSeq = cell(sizePicturePerSeq, sizeBlockType);

			% collect data in block loop.
			while ~feof(tmpFileAnalInfoPerSeq)
                tmpStrLine = fgetl(tmpFpClass);
                if tmpStrLine = 'P'
                	tmpIdxPic = str2num(tmpStrLine(end-1:end));
                elseif tmpStrLine ~= 'P' && tmpStrLine ~= 'b'
                	[tmpDepth, tmpMseEnc, tmpMseRef0, tmpMseRef1] = sscanf(tmpStrLine, '%d\t%d\t%d\t%d\n');
                	tmpAnalInfoPerSeq{tmpIdxPic+1, tmpDepth+1} = {tmpAnalInfoPerSeq{tmpIdxPic+1, tmpDepth+1}; [tmpMseEnc, tmpMseRef0, tmpMseRef1]};
                end
			end

			% collect memory.
			tmpAnalInfoPerClass{idxClass, idxQP} = tmpAnalInfoPerSeq;
		end
	end

	% collect memory.
	analInfo{idxClass} = tmpAnalInfoPerClass;
end