%=====================================
% Function: This script computes BD-rate gain of libvc under conditions
%           with different watching durations and/or start instants.
% File: computeLibvcGainInDiffDurationStart.m
% Author: ylonge from ZJU-MCL.
% Date: 2018-08-13.
% Update:
%=====================================

%%=====================================
% Definition Area.
% =====================================
% switches.
switchReadSeqList = 1;
switchReadFile = 1;
switchComputeGainLibvcWithDiffDuration = 1;
% constant values.
AVS_FORMAT = 1;
HEVC_FORMAT = 2;
ANCHOR_IDX = 1;
TEST_IDX = 2;
%fileFormat = AVS_FORMAT;
fileSeqList = '.\seqList_AVS3.txt';
pathBase = '.\copy\';
pathLogFile = [{'\anchor_avs2\'} {'\anchor_avs2_gb\'};{'\anchor_hm\'} {'\libvc\'}] ;
fileSubFix = '_enc.txt';
filePreFix = '\log_';
qpList = [27 32 38 45; 22 27 32 37];
sizeQpList = size(qpList, 2);

%%=====================================
% Read Sequence List Area.
% =====================================
if switchReadSeqList
% prepare values.
tmpFpSeqList = fopen(fileSeqList, 'r');
% prepare memory.
seqList = [];
intraPeriod = [];
% read file in loop.
while ~feof(tmpFpSeqList)
    tmpStrLine = fgetl(tmpFpSeqList);
    if tmpStrLine(1) == '#'
    	continue;
    end
    tmpStrSeq = sscanf(tmpStrLine, '%s\t%*[^#\n\r]');
    seqList = [seqList; {tmpStrSeq}];
    tmpIntraPeriod = sscanf(tmpStrLine, '%*s\t%*s\t%*f\t%*f\t%*f\t%*f\t%f%*[^#\n\r]');
    intraPeriod = [intraPeriod; tmpIntraPeriod];
end
% collect values.
sizeSeqList = length(seqList);

end % end of read sequence list.
%%=====================================
% Read File Area.
% =====================================
if switchReadFile
% prepare memory.
pictureDataAllPlatForm = cell(2, 1);
for fileFormat = 1:2
% prepare memory.
pictureData = cell(sizeSeqList, sizeQpList);
for idxSeq = 1: sizeSeqList
	for idxQP = 1: sizeQpList
		% prepare memory.
		tmpPictureDataPerSeqQP = cell(2, 2);
		for idxAnchorTest = 1:2
			% prepare values.
			tmpFileLog = [pathBase seqList{idxSeq} pathLogFile{fileFormat, idxAnchorTest} filePreFix seqList{idxSeq} '_' num2str(qpList(fileFormat, idxQP)) fileSubFix];
			tmpFpLog = fopen(tmpFileLog, 'r');
            if ferror(tmpFpLog)
                disp(['file ' tmpFileLog ' cannot be open\n']);
            end
			% prepare memory.
			tmpPictureDataPerType = [];
			tmpBasePictureDataPerType = [];

			% read data.
			if fileFormat == AVS_FORMAT
				% prepare values.
				tmpCountGBPic = 0;
				tmpIsNewGBPic = 0;
				while ~feof(tmpFpLog)
				    tmpStrLine = fgetl(tmpFpLog);
				    % check whether the line is empty.
				    if isempty(tmpStrLine)
				    	continue;
				    end
				    % check whether the background training is closed and a new GB picture is generated.
				    if tmpStrLine(1:7) == 'train c'
				    	tmpCountGBPic = tmpCountGBPic + 1;
				    	tmpIsNewGBPic = 1;
				    	continue;
				    end
				    % get the picture type.
				    tmpStrPicType = sscanf(tmpStrLine, '%*[^(](%[^)]%*[^#\n\r]');
                    if isempty(tmpStrPicType)
				    	continue;
				    end
				    tmpPicType = -1;
				    if tmpStrPicType(end) == 'I' || tmpStrPicType(end) == 'S'
				    	tmpPicType = 0;
				    elseif (tmpStrPicType(end) == 'B' && tmpStrPicType(1) ~= 'G') || tmpStrPicType(end) == 'P' || tmpStrPicType(end) == 'F'
				    	tmpPicType = 1;
				    elseif tmpStrPicType(end) == 'G'
				    	tmpPicType = 2;
				    	tmpCountGBPic = tmpCountGBPic + 1;
				    	tmpIsNewGBPic = 1;
				    elseif tmpStrPicType(end) == 'B' && tmpStrPicType(1) == 'G'
				    	tmpPicType = 3;
				    end
				    if tmpPicType == -1
				    	continue;
				    end
				    % get the poc.
				    tmpStrPoc = sscanf(tmpStrLine, '%[^(]%*[^#\n\r]');
				    if isempty(tmpStrPoc)
				    	continue;
				    end
				    if length(tmpStrPoc) < 6
				    	tmpPoc = str2num(tmpStrPoc);
				    else
				    	tmpPoc = str2num(tmpStrPoc(2:5));
				    end

				    tmpData  = sscanf(tmpStrLine, '%*[^)]) %f %*f %f %f %f %*[^#\n\r]');
				    if ~isempty(tmpData)
					    tmpBit = tmpData(1);
					    tmpPSNRY = tmpData(2);
					    tmpPSNRU = tmpData(3);
					    tmpPSNRV = tmpData(4);
					    if tmpPicType == 3 && tmpIsNewGBPic == 1 % GB
					    	tmpBasePictureDataPerType = [tmpBasePictureDataPerType; tmpCountGBPic tmpBit];
					    	tmpIsNewGBPic = 0;
                        elseif tmpPicType == 2 && tmpIsNewGBPic == 1% G
					    	tmpBasePictureDataPerType = [tmpBasePictureDataPerType; tmpCountGBPic tmpBit];
					    	tmpPictureDataPerType = [tmpPictureDataPerType; tmpPoc tmpCountGBPic 0 tmpPSNRY tmpPSNRU tmpPSNRV];
					    	tmpIsNewGBPic = 0;
                        elseif tmpPicType == 0 || tmpPicType == 1 % I S B P F
					    	tmpPictureDataPerType = [tmpPictureDataPerType; tmpPoc tmpCountGBPic tmpBit tmpPSNRY tmpPSNRU tmpPSNRV];
					    end
					end
				end
            elseif fileFormat == HEVC_FORMAT
				tmpIsLib = 0;
				tmpCountLibPic = 0;
				while ~feof(tmpFpLog)
				    tmpStrLine = fgetl(tmpFpLog);
                    % check whether the line is empty.
                    if isempty(tmpStrLine)
                        continue;
                    end
				    if tmpStrLine(1) == '#' && tmpIsLib == 0
				    	tmpIsLib = 1;
				    	continue;
				    elseif tmpStrLine(1) == '#' && tmpIsLib == 1
				    	tmpIsLib = 0;
				    	continue;
				    end
				    if tmpStrLine(1:3) ~= 'POC'
				    	continue;
				    end
				    if tmpIsLib == 0
				    	tmpStrPicType = sscanf(tmpStrLine, '%*[^(]( %[^,]%*[^\n\r]');
            			tmpData = sscanf(tmpStrLine, 'POC %f %*[^)]) %f bits [Y %f dB U %f dB V %f dB] %*[^]]] [L0 %f %*[^\n\r]');
            			if isempty(tmpData)
                            continue;
                        end
                        tmpPoc = tmpData(1);
	            		tmpBit = tmpData(2);
	            		tmpPSNRY = tmpData(3);
	            		tmpPSNRU = tmpData(4);
	            		tmpPSNRV = tmpData(5);
	            		if tmpStrPicType == 'P-SLICE'
	            			tmpRefLibIdx = tmpData(6);
	            		else
	            			tmpRefLibIdx = -1;
	            		end
					    tmpPictureDataPerType = [tmpPictureDataPerType; tmpPoc tmpRefLibIdx tmpBit tmpPSNRY tmpPSNRU tmpPSNRV];
            		else
            			tmpData = sscanf(tmpStrLine, 'POC %f %*[^)]) %f bits [Y %f dB U %f dB V %f %*[^\n\r]');
            			tmpPoc = tmpData(1);
	            		tmpBit = tmpData(2);
					    tmpBasePictureDataPerType = [tmpBasePictureDataPerType; tmpPoc tmpBit];
            		end
				end
			end % end read data.
			% collect data.
			tmpPictureDataPerSeqQP{idxAnchorTest, 1} = tmpBasePictureDataPerType;
			tmpPictureDataPerSeqQP{idxAnchorTest, 2} = tmpPictureDataPerType;
            % clear memory
            fclose(tmpFpLog);
		end % end of anchor test.
		% collect data.
		pictureData{idxSeq, idxQP} = tmpPictureDataPerSeqQP;
	end % end of QP.
end % end of Seq.
pictureDataAllPlatForm{fileFormat} = pictureData;
end % end of fileFormat.
				
end % end of read file.

%%=====================================
% Compute BD-rate gain under conditions with different durations.
% =====================================
if switchComputeGainLibvcWithDiffDuration
% prepare memory.
gainAllSeqAllPlatForm = cell(2, 1);
for fileFormat = 1:2
% prepare values.
durationList = 1:30; % seconds unit.
sizeDuration = length(durationList);
pictureData = pictureDataAllPlatForm{fileFormat};
% prepare memory.
gainAllSeq = cell(sizeSeqList, 1);
for idxSeq = 1: sizeSeqList
	% prepare values.
	tmpIntraPeriod = intraPeriod(idxSeq);
	% prepare memory.
	tmpGainPerSeq = zeros(sizeDuration, 3);
	for idxDuration = 1: sizeDuration
		% prepare values.
		tmpSizePic = size(pictureData{idxSeq, 1}{1, 2}, 1);
		if tmpIntraPeriod * idxDuration > tmpSizePic
			break;
		end
		% prepare memory.
		tmpBitPsnr = cell(2, 1); % 1 for anchor, 2 for test.
		for idxAnchorTest = 1: 2
			% prepare memory.
			tmpBitPsnrAnchorTest = zeros(sizeQpList, 4);
			% process data.
			for idxQP = 1: sizeQpList
                % prepare vales.
                tmpPocFirst = 0;
                tmpPocEnd = tmpPocFirst + tmpIntraPeriod * idxDuration - 1;
				% prepare memory.
				tmpAveragePictureData = [];
				while tmpPocEnd < tmpSizePic % check whether picturen in given duration is in the sequence range.
					% prepare picture data.
					tmpPictureDataPerSeqQPAnchorTest = pictureData{idxSeq, idxQP}{idxAnchorTest, 2};
					tmpBasePictureDataPerSeqQPAnchorTest = pictureData{idxSeq, idxQP}{idxAnchorTest, 1};
					% fetch picture data in range.
					tmpPictureDataInRange = tmpPictureDataPerSeqQPAnchorTest(tmpPictureDataPerSeqQPAnchorTest(:, 1)>=tmpPocFirst ...
						& tmpPictureDataPerSeqQPAnchorTest(:, 1)<=tmpPocEnd, :);
					% fetch base picture data in range.
					if isempty(tmpBasePictureDataPerSeqQPAnchorTest)
						tmpBasePictureBitInRange = 0;
					else
						tmpBasePicturePocInRange = unique(tmpPictureDataInRange(tmpPictureDataInRange(:, 2)>=0, 2));
                        tmpBasePictureIdxInRange = zeros(size(tmpBasePictureDataPerSeqQPAnchorTest, 1), 1);
                        for idx = 1: length(tmpBasePicturePocInRange)
                            tmpBasePictureIdxInRange = tmpBasePictureIdxInRange | tmpBasePictureDataPerSeqQPAnchorTest(:, 1) == tmpBasePicturePocInRange(idx);
                        end
						tmpBasePictureBitInRange = tmpBasePictureDataPerSeqQPAnchorTest(tmpBasePictureIdxInRange, 2);
					end
					% compute average data in range.
					tmpSumPictureDataInRange = sum(tmpPictureDataInRange(:, 3:end), 1);
					tmpSumPictureDataInRange(1) = tmpSumPictureDataInRange(1) + sum(tmpBasePictureBitInRange);
					tmpAveragePictureDataInRange = tmpSumPictureDataInRange ./ (tmpPocEnd - tmpPocFirst + 1);
					% collect data.
					tmpAveragePictureData = [tmpAveragePictureData; tmpAveragePictureDataInRange];
					% update range.
					tmpPocFirst = tmpPocFirst + tmpIntraPeriod;
					tmpPocEnd = tmpPocFirst + tmpIntraPeriod * idxDuration - 1;
				end
				% compute average data for given duration.
				tmpAverageBitPsnr = mean(tmpAveragePictureData, 1);
				tmpAverageBitPsnr(1) = tmpAverageBitPsnr(1) / 1000 / (1 / round(tmpIntraPeriod/10) / 10);
				tmpBitPsnrAnchorTest(idxQP, :) = tmpAverageBitPsnr;
			end % end of QP.
			% collect data.
			tmpBitPsnr{idxAnchorTest} = tmpBitPsnrAnchorTest;
		end
		% compute BD-rate gain and collect data.
		tmpGainPerSeq(idxDuration, :) = bdRateComparation(tmpBitPsnr{1}, tmpBitPsnr{2});
	end
	% collect data.
	gainAllSeq{idxSeq} = tmpGainPerSeq;
end
gainAllSeqAllPlatForm{fileFormat} = gainAllSeq;
end % end of fileFormat.
end % end of switchComputeGainLibvcWithDiffDuration.

%%=====================================
% Clear temporary and index values after all process is done.
% =====================================
% clear memory.
clear tmp* idx*;