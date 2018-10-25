%% This script is used to read library information for all sequences with all QPs and delta QPs.
flagReadSeqListTest = 1;
flagReadSeqListTrain = 0;
flagReadLibInfo = 0;
flagReadLogFile = 0;
flagPrepareContentInfo = 0;
flagFitForBestDqpEachLib = 0;
flagCheckResultValidity = 0;
flagModelVerify = 1;
flagPlotIdealCurve = 0;

% prepare basic information.
listQP = 22:5:37;
listDqp = 0:-1:-11;
numQP = length(listQP);
numDqp = length(listDqp);
numReso = 4; % for 1920-1280-720-416

%% read sequences list file of CTC.
if flagReadSeqListTest
	% display information.
	disp('reading sequence list of CTC.\n');

    tempRcVecSeqsB = readSeqList('.\\SeqList_B.txt');
    tempRcVecSeqsC = readSeqList('.\\SeqList_C.txt');
    tempRcVecSeqsD = readSeqList('.\\SeqList_D.txt');
    tempNumVecSeqsB = size(tempRcVecSeqsB, 1);
    tempNumVecSeqsC = size(tempRcVecSeqsC, 1);
    tempNumVecSeqsD = size(tempRcVecSeqsD, 1);

    % collect data.
    cellListSeq_Test = {tempRcVecSeqsB ; {NaN}; tempRcVecSeqsC; tempRcVecSeqsD};
    numSeq_Test = [tempNumVecSeqsB; 0; tempNumVecSeqsC; tempNumVecSeqsD];

    % set path.
    path = '.\log-ctc\';
    pathYuv = '\\mcl\MCL_Space\TestSeq\';
    
    % clear temp data.
    clear tempRcVecSeqsB tempRcVecSeqsC tempRcVecSeqsD tempNumVecSeqsB tempNumVecSeqsC tempNumVecSeqsD;
end

%% read sequences list file of train sequences.
if flagReadSeqListTrain
    % display information.
	disp('reading train sequence list.\n');

    tempRcVecSeqs = readSeqList('.\\seqlist-train-all-resolution.txt');
    tempNumVecSeqs = size(tempRcVecSeqs, 1);

    % collect data.
    cellListSeq_Train = {tempRcVecSeqs(1:11, :); tempRcVecSeqs(12:27, :); tempRcVecSeqs(28:32, :); tempRcVecSeqs(33:37, :)};
    numSeq_Train = [size(cellListSeq_Train{1},1); size(cellListSeq_Train{2},1); size(cellListSeq_Train{3},1); size(cellListSeq_Train{4},1)];
    
    % set path.
    path = '.\log-train\';
    pathYuv = '\\mcl\MCL_Space\TestSeq\';

    % clear temp data.
    clear tempRcVecSeqs tempNumVecSeqs;
end

%% read library information file.
if flagReadLibInfo
	if flagReadSeqListTrain
		cellListSeq = cellListSeq_Train;
		numAllSeq = numSeq_Train;
	else
		cellListSeq = cellListSeq_Test;
		numAllSeq = numSeq_Test;
	end
	% display information.
	disp('reading library information of sequences.\n');
	% collect data.
	cellLibInfo = cell(numReso, 1);
	cellEncAnalInfo = cell(numReso, 1);
	for idxReso = 1: numReso
		if numAllSeq(idxReso) < 1
			continue;
		end
		% prepare data.
		tempNumAllSeq = numAllSeq(idxReso);
		% collect data.
		tempLibInfoPerReso = cell(tempNumAllSeq, 6);
		tempEncAnalInfoPerReso = cell(tempNumAllSeq, 1);
		for idxSeq = 1: tempNumAllSeq
			% for one sequence, libInfo is same for all qp and dqp, except the area of inter prediction.
			tempFileLibInfo = [path 'analysis\libraryInfo_' cellListSeq{idxReso}{idxSeq, 2} '.txt'];
			tempLibInfo = readLibInfo( tempFileLibInfo );

			tempListOrgLib = cell2mat(tempLibInfo(4));
			tempListKeyPic = cell2mat(tempLibInfo(11));
			tempListRefLib = cell2mat(tempLibInfo(12));
			tempListLibFreq = cell2mat(tempLibInfo(5));
			tempListCenterCost = cell2mat(tempLibInfo(7));
			tempListClusterCost = cell2mat(tempLibInfo(8));

			tempLibInfoPerReso(idxSeq, 1) = {tempListCenterCost};
			tempLibInfoPerReso(idxSeq, 2) = {tempListClusterCost};
			tempLibInfoPerReso(idxSeq, 3) = {tempListLibFreq};
			tempLibInfoPerReso(idxSeq, 4) = {tempListKeyPic};
			tempLibInfoPerReso(idxSeq, 5) = {tempListOrgLib};
			tempLibInfoPerReso(idxSeq, 6) = {tempListRefLib};

			% read encode analysis information for all qp and dqp.
			% format-- poc	 InterPercent	skipPercent	 interPxlCount	 interResiVar	 interTrCoeffVar	 intraPxlCount	 intraResiVar	 intraTrCoeffVar
			tempEncAnalInfoPerQP = cell(numQP, 1);
            if flagReadSeqListTest
                for idxQP = 1: numQP
                    tempStrQP = num2str(listQP(idxQP));
                    tempEncAnalInfoPerDqp = cell(numDqp, 1);
                    for idxDqp = 1: numDqp
                        tempStrDqp = num2str(-listDqp(idxDqp));
                        tempFileEncAnalInfo = [path 'encAnalysis\' cellListSeq{idxReso}{idxSeq, 2} '\tempAnalysisInfo_' tempStrQP '_' tempStrDqp '.txt'];
                        tempEncAnalInfo = readEncAnalysisInfo( tempFileEncAnalInfo );
                        tempEncAnalInfoPerDqp{idxDqp} = tempEncAnalInfo;
                    end
                    % collect original residue.
                    tempFileEncAnalInfo = [path 'encAnalysis\' cellListSeq{idxReso}{idxSeq, 2} '\tempAnalysisInfo_' tempStrQP '_' tempStrQP '.txt'];
                    tempEncAnalInfo = readEncAnalysisInfo( tempFileEncAnalInfo );
                    tempEncAnalInfoPerDqp{numDqp+1} = tempEncAnalInfo;
                    tempEncAnalInfoPerQP{idxQP} = tempEncAnalInfoPerDqp;
                end
                % collect original redisue variance.
                tempFileEncAnalInfo = [path 'encAnalysis\' cellListSeq{idxReso}{idxSeq, 2} '\tempAnalysisInfo_0_0.txt'];
                tempEncAnalInfo = readEncAnalysisInfo( tempFileEncAnalInfo );
                tempEncAnalInfoPerQP{numQP + 1} = tempEncAnalInfo;
            end
			tempEncAnalInfoPerReso{idxSeq} = tempEncAnalInfoPerQP;
		end % seq
		% collect data.
		cellLibInfo{idxReso} = tempLibInfoPerReso;
		cellEncAnalInfo{idxReso} = tempEncAnalInfoPerReso;
	end

	% collect data for train or test.
	if flagReadSeqListTrain
		cellLibInfo_Train = cellLibInfo;
		cellEncAnalInfo_Train = cellEncAnalInfo;
	elseif flagReadSeqListTest
		cellLibInfo_Test = cellLibInfo;
		cellEncAnalInfo_Test = cellEncAnalInfo;
	end

	% clear temp data.
	clear temp* idx*;
end

%% read log file.
if flagReadLogFile
	subFlagReadSeqLog = 0;
	subFlagReadKeyLog = 0;
	subFlagReadBitPsnrEachLib = 1;

	%% read sequence encoding log file for bitrate and psnr.
	if subFlagReadSeqLog
		if flagReadSeqListTrain
			cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
		else
			cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
		end
		% collect data.
		gainYAllSeq = cell(numReso, 1);

		for idxReso = 1: numReso
			if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumAllSeq = numAllSeq(idxReso);
			% collect data.
			tempGainYAllSeq = zeros(tempNumAllSeq, numDqp);
			for idxSeq = 1: tempNumAllSeq
				% read anchor.
				tempBitPsnrAnchor = zeros(4, 5);
				for idxQP = 1: numQP
					tempStrQP = num2str(listQP(idxQP));
					tempFileAnchor = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_anchor_gop16_enc.txt'];
					tempBitPsnrAnchor(idxQP, :) = readSumBitPsnr(tempFileAnchor);
				end

				% read libvc.
				% collect data.
				for idxDqp = 1: numDqp
					if(listDqp(idxDqp) <= 0)
						tempStrDqp = num2str(-listDqp(idxDqp));
					else
						tempStrDqp = ['0' num2str(listDqp(idxDqp))];
					end
					tempBitPsnrLib = zeros(4, 5);
					for idxQP = 1: numQP
						tempStrQP = num2str(listQP(idxQP));
						tempFileLib = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_libvc_gop16_fixQP_' tempStrDqp '_enc.txt'];
						tempBitPsnrLib(idxQP, :) = readSumBitPsnr(tempFileLib);
					end
					tempGain = bdRateComparation( tempBitPsnrAnchor, tempBitPsnrLib );
					tempGainYAllSeq(idxSeq, idxDqp) = tempGain(1);
	            end
			end
			% collect data.
			gainYAllSeq{idxReso} = tempGainYAllSeq;
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			gainYUVAllSeq_Train = gainYAllSeq;
		elseif flagReadSeqListTest
			gainYUVAllSeq_Test = gainYAllSeq;
		end

		% clear temp data.
		clear temp* idx*;
	end

	%% read sequence encoding log file for BD-rate of key frames.
	if subFlagReadKeyLog
		if flagReadSeqListTrain
			cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
		else
			cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
		end
		% collect data.
		gainKeyYAllSeq = zeros(numReso, 1);
		for idxReso = 1: numReso
			if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumAllSeq = numAllSeq(idxReso);
			% collect data.
			tempGainKeyYAllSeq = zeros(tempNumAllSeq, numDqp);
			for idxSeq = 1: tempNumAllSeq
				% read anchor.
				tempBitPsnrAnchor = zeros(4, 5);
				for idxQP = 1: numQP
					tempStrQP = num2str(listQP(idxQP));
					tempFileAnchor = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_anchor_gop16_enc.txt'];
					tempBitPsnrAnchor(idxQP, :) = readKeyFrameSumBitPsnr(tempFileAnchor);
				end

				for idxDqp = 1: numDqp
					if(listDqp(idxDqp) <= 0)
						tempStrDqp = num2str(-listDqp(idxDqp));
					else
						tempStrDqp = ['0' num2str(listDqp(idxDqp))];
					end
					tempBitPsnrLib = zeros(4, 5);
					for idxQP = 1: numQP
						tempStrQP = num2str(listQP(idxQP));
						tempFileLib = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_libvc_gop16_fixQP_' tempStrDqp '_enc.txt'];
						tempBitPsnrLib(idxQP, :) = readKeyFrameSumBitPsnr(tempFileLib);
					end
					tempGain = bdRateComparation( tempBitPsnrAnchor, tempBitPsnrLib );
					tempGainKeyYAllSeq(idxSeq, idxDqp) = tempGain(1);
				end
			end
			% collect data.
			gainKeyYAllSeq(idxReso) = tempGainKeyYAllSeq;
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			gainKeyYAllSeq_Train = gainKeyYAllSeq;
		elseif flagReadSeqListTest
			gainKeyYAllSeq_Test = gainKeyYAllSeq;
		end

		% clear temp data.
		clear temp* idx*;
	end

	%% read average bitrate and psnr of frames for each library picture.
	if subFlagReadBitPsnrEachLib
		if flagReadSeqListTrain
			cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
		else
			cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
		end
		% display information.
		disp('reading bitrate and PSNR for each cluster of one library picture.\n');
		% collect data.
		gainYEachLibAllSeq = cell(numReso, 2); % 1 for sequence-wise, 2 for key-wise.
		scalePicDivKeyRDAllSeq = cell(numReso, 2);
		numFrameEachLibAllSeq = cell(numReso, 1);
		bitPsnrEachLibAllSeq = cell(numReso, 1);
		bitPsnrAllPicAllLibAllSeq = cell(numReso, 1);
		for idxReso = 1: numReso
			if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempGainYEachLibAllSeq = cell(tempNumSeq, numDqp);
			tempGainYEachLibAllSeqKeyFrames = cell(tempNumSeq, numDqp);
		    tempNumFrameEachLibAllSeq = cell(tempNumSeq, 1);
		    tempScalePicDivKeyRDAnchorAllSeq = cell(tempNumSeq, 1);
		    tempScalePicDivKeyRDLibvcAllSeq = cell(tempNumSeq, 1);
		    tempBitPsnrEachLibAllSeq = cell(tempNumSeq, 2); % 1 for anchor and 2 for libvc.
			tempBitPsnrAllPicAllLibAllSeq = cell(tempNumSeq, 3);
			for idxSeq = 1: tempNumSeq
		    	% prepare data.
		        tempFrameRate = cellListSeq{idxReso}{idxSeq, 3};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
		    	tempNumLibPic = length(tempListOrgLib);
		    	tempNumKeyPic = length(tempListKeyPic);
				% read anchor.
				% collect data.
				tempBitPsnrAnchorEachLib = cell(numQP, 1);
				tempBitPsnrAnchorAllLib = zeros(numQP* tempNumLibPic, 4);
				tempBitPsnrAnchorAllLibKeyFrames = zeros(numQP* tempNumLibPic, 4);
				tempScalePicDivKeyRDAnchorAllQPAllKey = zeros(numQP * tempNumKeyPic, 4);
                tempBitPsnrAllPicAllLibPerSeq = cell(numQP, 1);
				for idxQP = 1: numQP
					tempStrQP = num2str(listQP(idxQP));
					tempFileAnchor = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_anchor_gop16_enc.txt'];
					[tempBitPsnrAnchor, tempBitPsnrAnchorKeyFrames, tempScalePicDivKeyRD, tempListNumFrames, tempListBitDistAllPicAllLib] = readBitPsnrSingleLib(tempFileAnchor, tempListOrgLib, tempListKeyPic, tempListRefLib, tempFrameRate, 0);
					tempBitPsnrAnchorAllLib(idxQP: 4: end, :) = tempBitPsnrAnchor;
					tempBitPsnrAnchorAllLibKeyFrames(idxQP: 4: end, :) = tempBitPsnrAnchorKeyFrames;
					tempScalePicDivKeyRDAnchorAllQPAllKey(((idxQP-1)*tempNumKeyPic+1):idxQP*tempNumKeyPic, :) = tempScalePicDivKeyRD;
                    tempBitPsnrAllPicAllLibPerSeq{idxQP} = tempListBitDistAllPicAllLib;
					
					% collect bit and psnr for best gain compute based on the train model.
					tempBitPsnrAnchorEachLibPic = cell(tempNumLibPic, 2); % 1 for bit and psnr and 2 for frame number.
					for idxLibPic = 1: tempNumLibPic
						tempBitPsnrAnchorEachLibPic{idxLibPic, 1} = tempBitPsnrAnchor(idxLibPic, :);
						tempBitPsnrAnchorEachLibPic{idxLibPic, 2} = tempListNumFrames(idxLibPic);
					end
					tempBitPsnrAnchorEachLib{idxQP} = tempBitPsnrAnchorEachLibPic;
				end
				% collect data.
				tempNumFrameEachLibAllSeq(idxSeq) = {tempListNumFrames};
				tempBitPsnrEachLibAllSeq{idxSeq, 1} = tempBitPsnrAnchorEachLib;
                tempBitPsnrAllPicAllLibAllSeq{idxSeq, 1} = tempBitPsnrAllPicAllLibPerSeq;

				% read library data.
				% collect data.
				tempScalePicDivKeyRDLibvcAllQPAllKey = cell(numQP, 1);
				tempBitPsnrLibvcEachLib = cell(numQP, 1);
				tempBitPsnrAllPicAllLibPerSeq = cell(numQP, 1);
                tempDiffBitPsnrAllPicAllLibPerSeq = cell(numQP, 1);
				for idxQP = 1: numQP
					tempScalePicDivKeyRDLibvcAllQPAllKey{idxQP} = cell(numDqp, 1);
					tempBitPsnrLibvcEachLib{idxQP} = cell(numDqp, 1);
					tempBitPsnrAllPicAllLibPerSeq{idxQP} = cell(numDqp, 1);
                    tempDiffBitPsnrAllPicAllLibPerSeq{idxQP} = cell(numDqp, 1);
				end
				for idxDqp = 1: numDqp
					% prepare path.
					if(listDqp(idxDqp) <= 0)
						tempStrDqp = num2str(-listDqp(idxDqp));
					else
						tempStrDqp = ['0' num2str(listDqp(idxDqp))];
					end
					% collect data.
					tempBitPsnrLibAllLib = zeros(numQP* tempNumLibPic, 4);
					tempBitPsnrLibAllLibKeyFrames = zeros(numQP* tempNumLibPic, 4);
					for idxQP = 1: numQP
						tempStrQP = num2str(listQP(idxQP));
						tempFileLib = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_libvc_gop16_fixQP_' tempStrDqp '_enc.txt'];
						[tempBitPsnrLib, tempBitPsnrLibKeyFrames, tempScalePicDivKeyRD, tempListNumFrames, tempListBitDistAllPicAllLib] = readBitPsnrSingleLib(tempFileLib, tempListOrgLib, tempListKeyPic, tempListRefLib, tempFrameRate, 1);
						tempBitPsnrLibAllLib(idxQP: 4: end, :) = tempBitPsnrLib;
						tempBitPsnrLibAllLibKeyFrames(idxQP: 4: end, :) = tempBitPsnrLibKeyFrames;
						tempScalePicDivKeyRDLibvcAllQPAllKey{idxQP}{idxDqp} = tempScalePicDivKeyRD;
						tempBitPsnrAllPicAllLibPerSeq{idxQP}{idxDqp} = tempListBitDistAllPicAllLib;
%                         for i = 1: length(tempListBitDistAllPicAllLib)
%                             tempListBitDistAllPicAllLib{i} = tempListBitDistAllPicAllLib{i} - tempBitPsnrAllPicAllLibAllSeq{idxSeq, 1}{idxQP}{i};
%                         end
%                         tempDiffBitPsnrAllPicAllLibPerSeq{idxQP}{idxDqp} = tempListBitDistAllPicAllLib;

						% collect bit and psnr.
						tempBitPsnrLibvcEachLibPic = cell(tempNumLibPic, 2); % 1 for bit and psnr and 2 for frame number.
						for idxLibPic = 1: tempNumLibPic
							tempBitPsnrLibvcEachLibPic{idxLibPic, 1} = tempBitPsnrLib(idxLibPic, :);
							tempBitPsnrLibvcEachLibPic{idxLibPic, 2} = tempListNumFrames(idxLibPic);
						end
						tempBitPsnrLibvcEachLib{idxQP}{idxDqp} = tempBitPsnrLibvcEachLibPic;
					end

					% compute BD-rate for each group of library picture.
					tempGain = zeros(tempNumLibPic, 3);
					tempGainKeyFrames = zeros(tempNumLibPic, 3);
					for idxLibPic = 1: tempNumLibPic
						% BD-rate for all frames.
						tempBitPsnrAnchor = tempBitPsnrAnchorAllLib(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
						tempBitPsnrLib = tempBitPsnrLibAllLib(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
						tempGain(idxLibPic, :) = bdRateComparation( tempBitPsnrAnchor, tempBitPsnrLib );
						% BD-rate for key frames.
						tempBitPsnrAnchorKeyFrames = tempBitPsnrAnchorAllLibKeyFrames(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
						tempBitPsnrLibKeyFrames = tempBitPsnrLibAllLibKeyFrames(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
						tempGainKeyFrames(idxLibPic, :) = bdRateComparation( tempBitPsnrAnchorKeyFrames, tempBitPsnrLibKeyFrames );
					end
					tempGainYEachLibAllSeq(idxSeq, idxDqp) = {tempGain(:, 1)};
					tempGainYEachLibAllSeqKeyFrames(idxSeq, idxDqp) = {tempGainKeyFrames(:, 1)};
				end
				% collect data.
				tempScalePicDivKeyRDAnchorAllSeq{idxSeq, 1} = tempScalePicDivKeyRDAnchorAllQPAllKey;
				tempScalePicDivKeyRDLibvcAllSeq{idxSeq, 2} = tempScalePicDivKeyRDLibvcAllQPAllKey;
				tempBitPsnrEachLibAllSeq{idxSeq, 2} = tempBitPsnrLibvcEachLib;
				tempBitPsnrAllPicAllLibAllSeq{idxSeq, 2} = tempBitPsnrAllPicAllLibPerSeq;
                tempBitPsnrAllPicAllLibAllSeq{idxSeq, 3} = tempDiffBitPsnrAllPicAllLibPerSeq;
			end
			gainYEachLibAllSeq{idxReso, 1} = tempGainYEachLibAllSeq;
			gainYEachLibAllSeq{idxReso, 2} = tempGainYEachLibAllSeqKeyFrames;
			scalePicDivKeyRDAllSeq{idxReso, 1} = tempScalePicDivKeyRDAnchorAllSeq;
			scalePicDivKeyRDAllSeq{idxReso, 2} = tempScalePicDivKeyRDLibvcAllSeq;
			numFrameEachLibAllSeq{idxReso} = tempNumFrameEachLibAllSeq;
			bitPsnrEachLibAllSeq{idxReso} = tempBitPsnrEachLibAllSeq;
			bitPsnrAllPicAllLibAllSeq{idxReso} = tempBitPsnrAllPicAllLibAllSeq;
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			gainYEachLibAllSeq_Train = gainYEachLibAllSeq;
			numFrameEachLibAllSeq_Train = numFrameEachLibAllSeq;
			scalePicDivKeyRDAllSeq_Train = scalePicDivKeyRDAllSeq;
			bitPsnrEachLibAllSeq_Train = bitPsnrEachLibAllSeq;
			bitPsnrAllPicAllLibAllSeq_Train = bitPsnrAllPicAllLibAllSeq;
		elseif flagReadSeqListTest
			gainYEachLibAllSeq_Test = gainYEachLibAllSeq;
			numFrameEachLibAllSeq_Test = numFrameEachLibAllSeq;
			scalePicDivKeyRDAllSeq_Test = scalePicDivKeyRDAllSeq;
			bitPsnrEachLibAllSeq_Test = bitPsnrEachLibAllSeq;
			bitPsnrAllPicAllLibAllSeq_Test = bitPsnrAllPicAllLibAllSeq;
		end

		% clear temp data.
		clear temp* idx*;
	end
end

if flagPrepareContentInfo
	subFlagComputeContentInfo = 0;
	subFlagComputeDiffPic = 1;

    %% read sequence yuv to compute content information of picutures.
	if subFlagComputeContentInfo
		% prepare information.
		if flagReadSeqListTrain
			cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
		else
			cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
		end
		tempBlockSize = 16;
        tempSearchRange = 32;
		% display information.
		disp('computing content information from sequences.\n');

		% collect data.
		oriContentInfoAllSeq = cell(numReso, 1);
		diffPicAllSeq = cell(numReso, 1);
		for idxReso = 1: numReso
			if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempOriContentInfoAllSeq = cell(tempNumSeq, 2);
			tempDiffPicAllSeq = cell(tempNumSeq, 1);
	        
			for idxSeq = 1: tempNumSeq
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumKeyPic = length(tempListKeyPic);
				tempNumLibPic = length(tempListOrgLib);

				tempFileYuv = [pathYuv cellListSeq{idxReso}{idxSeq, 1} '\' cellListSeq{idxReso}{idxSeq, 2} '.yuv'];
	            disp(tempFileYuv);

	            % collect data.
	            tempEstInfoAllKeyPic = cell(tempNumKeyPic, 3);
	            tempEstInfoAllLibPic = cell(tempNumKeyPic, 1);
	            tempDiffPicAllKeyPic = cell(tempNumKeyPic, 1);

				for idxKeyPic = 1: tempNumKeyPic
					% prepare data.
					idxLibPic = tempListRefLib(idxKeyPic) + 1;
	                tempOrgPocLibPic = tempListOrgLib(idxLibPic);
	                tempOrgPocIntraPic = tempListKeyPic(idxKeyPic);
	                tempYKeyPic = readYuv(tempFileYuv, tempOrgPocIntraPic);
	                if tempOrgPocLibPic == tempOrgPocIntraPic
	                    tempYLibPic = tempYKeyPic;
	                else
	                    tempYLibPic = readYuv(tempFileYuv, tempListOrgLib(idxLibPic));
	                end

	                % method2: average block difference based on motion compensation.
                    if tempOrgPocLibPic ~= tempOrgPocIntraPic
                        [ tempErrInterAver, tempErrInterTotal, tempErrIntraAver, tempErrIntraTotal, tempPercentIntra, ~, tempErrMatInter, tempErrMatIntra] = motionCompensation( tempYLibPic, tempYKeyPic, tempBlockSize, tempSearchRange );
                        tempEstResiEnergyInterAllKeyPic = {tempErrInterAver; tempErrInterTotal; tempErrMatInter};
                        tempEstResiEnergyIntraAllKeyPic = {tempErrIntraAver; tempErrIntraTotal; tempErrMatIntra};
                        tempPercentInterAllKeyPic = 1 - tempPercentIntra;
                        tempEstInfoAllKeyPic(idxKeyPic, :) = {tempEstResiEnergyInterAllKeyPic, tempEstResiEnergyIntraAllKeyPic, tempPercentInterAllKeyPic};
                    else
                        [ tempErrIntraAver, tempErrIntraTotal, tempErrMatIntra ] = intraDifference( tempYLibPic, tempBlockSize );
                        tempEstResiEnergyAllLibPic = {tempErrIntraAver; tempErrIntraTotal; tempErrMatIntra};
                        tempEstInfoAllLibPic(idxKeyPic) = {tempEstResiEnergyAllLibPic};
                    end
	            end
	            % collect data.
				tempOriContentInfoAllSeq(idxSeq, 1) = {tempEstInfoAllKeyPic};
				tempOriContentInfoAllSeq(idxSeq, 2) = {tempEstInfoAllLibPic};
			end
			% collect data.
			oriContentInfoAllSeq{idxReso} =tempOriContentInfoAllSeq;
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			oriContentInfoAllSeq_Train = oriContentInfoAllSeq;
		elseif flagReadSeqListTest
			oriContentInfoAllSeq_Test = oriContentInfoAllSeq;
		end

		% clear temp data.
		clear temp* idx*;
    end

    %% read sequence yuv to compute content information of picutures.
	if subFlagComputeDiffPic
		% prepare information.
		if flagReadSeqListTrain
			cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
		else
			cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
		end

		% collect data.
		diffPicAllSeq = cell(numReso, 1);
		for idxReso = 1: numReso
			if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempDiffPicAllSeq = cell(tempNumSeq, 1);
	        
			for idxSeq = 1: tempNumSeq
				tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempSeqLength = cellListSeq{idxReso}{idxSeq, 4};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumKeyPic = length(tempListKeyPic);
				tempNumLibPic = length(tempListOrgLib);

				tempFileYuv = [pathYuv cellListSeq{idxReso}{idxSeq, 1} '\' cellListSeq{idxReso}{idxSeq, 2} '.yuv'];
	            disp(tempFileYuv);

	            % collect data.
	            tempDiffPicAllKeyPic = cell(tempNumLibPic, 1);

	            for idxLibPic = 1: tempNumLibPic
	            	tempNumKeyPicPerLib = sum(tempListRefLib == (idxLibPic - 1));
	            	tempDiffPicPerLib = zeros(tempNumKeyPicPerLib, 1);
	            	tempCountKeyPicPerLib = 0;
					for idxKeyPic = 1: tempNumKeyPic
						if tempListRefLib(idxKeyPic) ~= idxLibPic - 1
							continue;
						end
						% prepare data.
		                tempOrgPocIntraPicCurr = tempListKeyPic(idxKeyPic);
		                if idxKeyPic == tempNumKeyPic
		                	tempOrgPocIntraPicAfter = tempSeqLength - 1;
		                else
		                	tempOrgPocIntraPicAfter = tempListKeyPic(idxKeyPic + 1) - 1;
		                end

		                % collect data.
		                tempAverageDiffPic = 0;
		                % compute difference.
		               	tempYPicAfter = readYuv(tempFileYuv, tempOrgPocIntraPicCurr);
		                for idxPoc = tempOrgPocIntraPicCurr : tempOrgPocIntraPicAfter - 1
		                	tempYPicCur = tempYPicAfter;
		                	tempYPicAfter = readYuv(tempFileYuv, idxPoc + 1);
		                	tempDiffPic = sum(abs(tempYPicAfter(:) - tempYPicCur(:))) / tempWidth / tempHeight;
		                	tempAverageDiffPic = tempAverageDiffPic + tempDiffPic;
		                end
		                tempAverageDiffPic = tempAverageDiffPic / (tempOrgPocIntraPicAfter - tempOrgPocIntraPicCurr);
		                tempCountKeyPicPerLib = tempCountKeyPicPerLib + 1;
		                tempDiffPicPerLib(tempCountKeyPicPerLib) = tempAverageDiffPic;
		            end
		            tempDiffPicAllKeyPic{idxLibPic} = tempDiffPicPerLib;
		        end
	            % collect data.
	            tempDiffPicAllSeq{idxSeq} = tempDiffPicAllKeyPic;
			end
			% collect data.
			diffPicAllSeq{idxReso} =tempDiffPicAllSeq;
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			diffPicAllSeq_Train = diffPicAllSeq;
		elseif flagReadSeqListTest
			diffPicAllSeq_Test = diffPicAllSeq;
		end

		% clear temp data.
		clear temp* idx*;
    end
end

%% check the result.
if flagCheckResultValidity
	subFlagFindBestDqpForAllSeqAllLibPicAllQp = 1;
	subFlagShowResAtAllQPCombination = 1;

	%% find the best delta QP for all sequences, all library pictures, all QP.
    if subFlagFindBestDqpForAllSeqAllLibPicAllQp
    	% display information.
		disp('finding the best delta QP for all sequences, all library pictures, all QP.\n');

		% prepare information.
		if flagReadSeqListTrain
			cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
		else
			cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
		end
		% collect data.
		bestDqpAndGainAllLibPicAllQp_BDrate = cell(numReso, 3); % 1 for sequence-wise, 2 for key-pic-wise. 3 for all bdrate at all prossible QP combinations.
		for idxReso = 1: numReso
			if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempBestDqpAllSeqAllLibPicAllQp_BDrate = [];
			tempBestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate = [];
	        tempBestGainAllSeqAllLibPic_BDrate = [];
	        tempBestGainAllSeqAllLibPicKeyFrames_BDrate = [];
	        tempGainAllQPCombinationAllSeq = cell(tempNumSeq, 1);

			for idxSeq = 1: tempNumSeq
		    	% prepare data.
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
		    	tempNumLibPic = length(tempListOrgLib);
		        tempFrameRate = cellListSeq{idxReso}{idxSeq, 3};
				% read anchor.
				% collect data.
				tempBitPsnrAnchorAllLib = cell(numQP, tempNumLibPic);
				tempBitPsnrAnchorAllLibKeyFrames = cell(numQP, tempNumLibPic);
				for idxQP = 1: numQP
					tempStrQP = num2str(listQP(idxQP));
					tempFileAnchor = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_anchor_gop16_enc.txt'];
					[tempBitPsnrAnchor, tempBitPsnrAnchorKeyFrames, tempScalePicDivKeyRD, tempListNumFrames, ~] = readBitPsnrSingleLib(tempFileAnchor, tempListOrgLib, tempListKeyPic, tempListRefLib, tempFrameRate, 0);
					for idxLibPic = 1: tempNumLibPic
						tempBitPsnrAnchorAllLib(idxQP, idxLibPic) = {tempBitPsnrAnchor(idxLibPic, 1: 4)};
						tempBitPsnrAnchorAllLibKeyFrames(idxQP, idxLibPic) = {tempBitPsnrAnchorKeyFrames(idxLibPic, 1: 4)};
					end
				end

				% collect data.
				tempBitPsnrLibPicAllLibAllDqp = cell(numQP, numDqp, tempNumLibPic);
				tempBitPsnrLibPicAllLibAllDqpKeyFrames = cell(numQP, numDqp, tempNumLibPic);
				for idxDqp = 1: numDqp
					% prepare path.
					if(listDqp(idxDqp) <= 0)
						tempStrDqp = num2str(-listDqp(idxDqp));
					else
						tempStrDqp = ['0' num2str(listDqp(idxDqp))];
					end
					% collect data.
					for idxQP = 1: numQP
						tempStrQP = num2str(listQP(idxQP));
						tempFileLib = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_libvc_gop16_fixQP_' tempStrDqp '_enc.txt'];
						[tempBitPsnrLib, tempBitPsnrLibKeyFrames, tempScalePicDivKeyRD, tempListNumFrames, ~] = readBitPsnrSingleLib(tempFileLib, tempListOrgLib, tempListKeyPic, tempListRefLib, tempFrameRate, 1);
						for idxLibPic = 1: tempNumLibPic
							tempBitPsnrLibPicAllLibAllDqp(idxQP, idxDqp, idxLibPic) = {tempBitPsnrLib(idxLibPic, 1: 4)};
							tempBitPsnrLibPicAllLibAllDqpKeyFrames(idxQP, idxDqp, idxLibPic) = {tempBitPsnrLibKeyFrames(idxLibPic, 1: 4)};
						end
					end
				end

				% collect data.
				tempGainAllQPCombinationPerSeq = cell(tempNumLibPic, 2); % 1 for gain and 2 for all bitrate/psnr.
				% compute BD-rate and find the dest delta QP.
				for idxLibPic = 1: tempNumLibPic
					% collect data.
					tempGainAllQPCombinationPerLib = zeros(numDqp*numDqp*numDqp*numDqp, 7); % 1-4 for dqp at four qp and 5-7 for gain of YUV.
					tempBitPsnrAllQPCombinationPerLib = cell(numDqp*numDqp*numDqp*numDqp, 1); % 1-4 for dqp at four qp and 5-7 for gain of YUV.
					% BD-rate for all frames.
					tempBestGain = Inf;
					tempBestGainKeyFrames = Inf;
					tempBestDqpIdx = zeros(numQP, 1);
					tempBestDqpIdxKeyFrames = zeros(numQP, 1);
					% temporal memory for QP 22, 27, 32, 37.
					tempBitPsnrAnchor = [tempBitPsnrAnchorAllLib{1, idxLibPic}; tempBitPsnrAnchorAllLib{2, idxLibPic}; tempBitPsnrAnchorAllLib{3, idxLibPic}; tempBitPsnrAnchorAllLib{4, idxLibPic};];
					tempBitPsnrAnchorKeyFrames = [tempBitPsnrAnchorAllLibKeyFrames{1, idxLibPic}; tempBitPsnrAnchorAllLibKeyFrames{2, idxLibPic}; tempBitPsnrAnchorAllLibKeyFrames{3, idxLibPic}; tempBitPsnrAnchorAllLibKeyFrames{4, idxLibPic};];
					tempBitPsnrLib = zeros(numQP, 4);
					tempBitPsnrLibKeyFrames = zeros(numQP, 4);

					for idxDqp1 = 1: numDqp
						tempBitPsnrLib(1, :) = tempBitPsnrLibPicAllLibAllDqp{1, idxDqp1, idxLibPic};
						tempBitPsnrLibKeyFrames(1, :) = tempBitPsnrLibPicAllLibAllDqpKeyFrames{1, idxDqp1, idxLibPic};
						for idxDqp2 = 1: numDqp
							tempBitPsnrLib(2, :) = tempBitPsnrLibPicAllLibAllDqp{2, idxDqp2, idxLibPic};
							tempBitPsnrLibKeyFrames(2, :) = tempBitPsnrLibPicAllLibAllDqpKeyFrames{2, idxDqp2, idxLibPic};
							for idxDqp3 = 1: numDqp
								tempBitPsnrLib(3, :) = tempBitPsnrLibPicAllLibAllDqp{3, idxDqp3, idxLibPic};
								tempBitPsnrLibKeyFrames(3, :) = tempBitPsnrLibPicAllLibAllDqpKeyFrames{3, idxDqp3, idxLibPic};
								for idxDqp4 = 1: numDqp
									tempBitPsnrLib(4, :) = tempBitPsnrLibPicAllLibAllDqp{4, idxDqp4, idxLibPic};
									tempBitPsnrLibKeyFrames(4, :) = tempBitPsnrLibPicAllLibAllDqpKeyFrames{4, idxDqp4, idxLibPic};

									% compute BD-rate.
									tempGain = bdRateComparation( tempBitPsnrAnchor, tempBitPsnrLib );
									if tempGain(1) < tempBestGain
										tempBestDqpIdx = [idxDqp1; idxDqp2; idxDqp3; idxDqp4];
	                                    tempBestGain = tempGain(1);
									end
									tempGainKeyFrames = bdRateComparation( tempBitPsnrAnchorKeyFrames, tempBitPsnrLibKeyFrames );
									if tempGainKeyFrames < tempBestGainKeyFrames
										tempBestDqpIdxKeyFrames = [idxDqp1; idxDqp2; idxDqp3; idxDqp4];
	                                    tempBestGainKeyFrames = tempGainKeyFrames(1);
									end
									% collect data.
									tempGainAllQPCombinationPerLib((((idxDqp1-1)*numDqp+idxDqp2-1)*numDqp+idxDqp3-1)*numDqp+idxDqp4, :) = [listDqp(idxDqp1) listDqp(idxDqp2) listDqp(idxDqp3) listDqp(idxDqp4) tempGain];
									tempBitPsnrAllQPCombinationPerLib{(((idxDqp1-1)*numDqp+idxDqp2-1)*numDqp+idxDqp3-1)*numDqp+idxDqp4} = tempBitPsnrLib;
								end
							end
						end
					end
					% collect data.
					tempBestDqpAllSeqAllLibPicAllQp_BDrate = [tempBestDqpAllSeqAllLibPicAllQp_BDrate; listDqp(tempBestDqpIdx)'];
					tempBestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate = [tempBestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate; listDqp(tempBestDqpIdxKeyFrames)'];
					tempBestGainAllSeqAllLibPic_BDrate = [tempBestGainAllSeqAllLibPic_BDrate; tempBestGain];
					tempBestGainAllSeqAllLibPicKeyFrames_BDrate = [tempBestGainAllSeqAllLibPicKeyFrames_BDrate; tempBestGainKeyFrames];
					tempGainAllQPCombinationPerSeq{idxLibPic, 1} = tempGainAllQPCombinationPerLib;
					tempGainAllQPCombinationPerSeq{idxLibPic, 2} = tempBitPsnrAllQPCombinationPerLib;
				end
				tempGainAllQPCombinationAllSeq{idxSeq} = tempGainAllQPCombinationPerSeq;
			end
			% collect data.
			bestDqpAndGainAllLibPicAllQp_BDrate{idxReso, 1} = [{tempBestDqpAllSeqAllLibPicAllQp_BDrate} {tempBestGainAllSeqAllLibPic_BDrate}];
			bestDqpAndGainAllLibPicAllQp_BDrate{idxReso, 2} = [{tempBestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate} {tempBestGainAllSeqAllLibPicKeyFrames_BDrate}];
			bestDqpAndGainAllLibPicAllQp_BDrate{idxReso, 3} = tempGainAllQPCombinationAllSeq;
		end

		if flagReadSeqListTrain
			bestDqpAndGainAllLibPicAllQp_BDrate_Train = bestDqpAndGainAllLibPicAllQp_BDrate;
		elseif flagReadSeqListTest
			bestDqpAndGainAllLibPicAllQp_BDrate_Test = bestDqpAndGainAllLibPicAllQp_BDrate;
		end

		% clear temp data.
		clear temp* idx*;
    end
end

%% This part 
if flagModelVerify
	subFlagReadVerifyData = 0;
	
	subFlagVerifyModelLib = 0;
	subFlagVerifyReferenceModel = 0;
	
	subFlagTrainModelCoeffLib = 0;
	subFlagTrainModelCoeffKey = 0;
	subFlagTrainScalePicDivKeyRD = 0;
	subFlagTestDQPModelANDComputeBDRate = 1;
	subFlagTestTrainModelANDComputeBDRate = 0;
	subFlagRdcostOnActualBitDistANDComputeBDRate = 0;
    subFlagTestTrainModelOnActualEncInfoANDComputeBDRate = 0;

	subFlagPlotActualRateDist = 0;
	subFlagPlotTheoryRateDistTrend = 0;
	subFlagCheckEncInfoVSEstiInfo = 0;

	%% read verify data of rate and psnr of library pictures and key pictures.
	if subFlagReadVerifyData
		% prepare information.
		if flagReadSeqListTrain
			cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
		else
			cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
		end
		% collect data.
		bitPsnrLibKeyAllSeq = cell(numReso, 1);
		for idxReso = 1: numReso
			if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempBitPsnrLibKeyAllSeq = cell(tempNumSeq, 1);
			% loop for all sequences.
			for idxSeq = 1: tempNumSeq
				% prepare data.
				tempListLibFreq = cellLibInfo{idxReso}{idxSeq, 3};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
				% collect data.
				tempBitPsnrLibKeyAllQp = cell(numQP, 1);
				% loop for each base QP.
				for idxQP = 1: numQP
					% prepare log path.
					tempStrQP = num2str(listQP(idxQP));
					% for each dqp.
					% collect data.
					tempBitPsnrALLDqp = cell(numDqp, 1);

					for idxDqp = 1: numDqp
						% collect data.
						tempBitPsnrAllLib = cell(tempNumLibPic, 2);

						% prepare log path.
						if(listDqp(idxDqp) <= 0)
							tempStrDqp = num2str(-listDqp(idxDqp));
						else
							tempStrDqp = ['0' num2str(listDqp(idxDqp))];
						end
						tempFileLib = [path cellListSeq{idxReso}{idxSeq, 2} '\Sequence_' tempStrQP '_libvc_gop16_fixQP_' tempStrDqp '_enc.txt'];

						for idxLibPic = 1: tempNumLibPic
							% prepare data.
							tempNumKeyPicPerLib = tempListLibFreq(idxLibPic);
							tempPocLibPic = idxLibPic - 1;

							% collect library picture.
							% datatype: dataPoc dataType dataQP dataBit dataPSNRY dataPSNRU dataPSNRV.
							tempBitPsnrAllLib(idxLibPic, 1) = {readBitPsnrWithPoc( tempFileLib, tempPocLibPic, 1 )};
							% collect key picture.
							tempListIntraPocPerLib = zeros(tempNumKeyPicPerLib, 1);
							tempCount = 1;
							for idxKeyPic = 1: tempNumKeyPic
								if tempListRefLib(idxKeyPic) ~= tempPocLibPic
									continue;
								else
									tempListIntraPocPerLib(tempCount) = tempListKeyPic(idxKeyPic);
									tempCount = tempCount + 1;
								end
							end
							tempBitPsnrAllLib(idxLibPic, 2) = {readBitPsnrWithPoc( tempFileLib, tempListIntraPocPerLib, 0)};
						end
						% collect data.
						tempBitPsnrALLDqp(idxDqp) = {tempBitPsnrAllLib};
					end
					% collect data.
					tempBitPsnrLibKeyAllQp(idxQP) = {tempBitPsnrALLDqp};
				end
				% collect data.
				tempBitPsnrLibKeyAllSeq(idxSeq) = {tempBitPsnrLibKeyAllQp};
	        end
	        % collect data.
	        bitPsnrLibKeyAllSeq{idxReso} = tempBitPsnrLibKeyAllSeq;
	    end

        if flagReadSeqListTrain
            bitPsnrLibKeyAllSeq_Train = bitPsnrLibKeyAllSeq;
        elseif flagReadSeqListTest
            bitPsnrLibKeyAllSeq_Test = bitPsnrLibKeyAllSeq;
        end

		% clear temp data.
		clear temp* idx*;
	end

	%% basic parameters.
	dIntra = 2/3;
	dInter = 5/6;
	dNorm = 0.5;

	%% distortion model.
	% general integral model.
	equationIntegral = @(a,b,x) (-0.5.*exp(-a.*x).*((a.*x+1).^2+1)+1+1./(1-exp(-x)).*(-0.5*exp(-(1+a).*x)).*(((1-b).*x+1).^2+1-exp(x).*((b.*x-1).^2+1)));
	equationEntropy = @(a,x) (exp(-a*x).*(log2((1-exp(-a*x))./(1-exp(-x)))+1+x.*(a*(1-exp(-x))+exp(-x))./(1-exp(-x))/log(2))-log2(1-exp(-a*x)));
	% simplified integral model.
	eqT = @(d,t) 1+t.*exp(-d.*t)./2./(1-exp(-t)).*(t-2*d*t-2);
	eqTLow = @(d,t) 1+t.*(t-2*d.*t-2)./2./exp(d.*t);

	%% verify the rate and distortion model of library pictures based on estimated residue energy.
	if subFlagVerifyModelLib
		% get basic data.
		if flagReadSeqListTrain
            cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Train;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Train;
        elseif flagReadSeqListTest
            cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
        end
        % collect data.
        fitResLibAllSeq = cell(numReso, 1);
        for idxReso	= 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempFitResLibAllSeq = cell(tempNumSeq, 1);
			for idxSeq = 1: tempNumSeq
				% prepare data.
	            tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempListLibFreq = cellLibInfo{idxReso}{idxSeq, 3};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
				tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
				tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length', tempEstInfoAllLibPic) > 0);

				% collect data. There are 3 models to test.
				tempFitResAllLibPic = zeros(tempNumLibPic, 8);

				% collect library picture data.
				for idxLibPic = 1: tempNumLibPic
					% prepare content information.
					tempEstResiEnergyLib = tempEstInfoAllLibPic{idxLibPic, 1}{1};
					% collect data.
					tempDistAllQP = zeros(numQP * numDqp, 1);
					tempRateAllQP = zeros(numQP * numDqp, 1);
					tempQstepAllQP = zeros(numQP * numDqp, 1);

					for idxQP = 1: numQP
						for idxDqp = 1: numDqp
							tempQstepAllQP((idxQP-1)*numDqp + idxDqp, 1) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(3);
							tempRateAllQP((idxQP-1)*numDqp + idxDqp, 1) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(4);
							tempDistAllQP((idxQP-1)*numDqp + idxDqp, 1) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(5);
						end
	                end
					% remove repeat data.
					[tempQstepAllQP, idxUnique, ~] = unique(tempQstepAllQP, 'rows');
					tempRateAllQP = tempRateAllQP(idxUnique) / tempWidth / tempHeight;
					tempDistAllQP = psnr2mse(tempDistAllQP(idxUnique));
					tempQstepAllQP = qp2qstep(tempQstepAllQP);
					% fit data.
					tempFlagShow = 0;
					if tempFlagShow
						subplot(2,1,1);
	                end
	                % distortion model.
	                eqDLibSimpBaseT = @(a,b,m,x) a.*eqT(dNorm,m.*x)+b;
					equationDLib = eqDLibSimpBaseT;
	                tempOpts = fitoptions(equationDLib);
	                tempOpts.StartPoint = [1 1 1];
	                tempOpts.Lower = [0 -Inf 0];
					[ tempFitRes, tempGof ] = computeEquationFit( tempQstepAllQP, tempDistAllQP, equationDLib, tempOpts, tempFlagShow );
					tempFitResAllLibPic(idxLibPic, 1:4) = [tempGof.rsquare tempFitRes.a tempFitRes.b tempFitRes.m];

	                if tempFlagShow
	                	subplot(2,1,2);
	                end
					% rate model.
					eqRLib = @(u,b,m,x) u.*equationEntropy(dNorm,m.*x)+b;
	                equationRLib = eqRLib;
	                tempOpts = fitoptions(equationRLib);
	                tempOpts.StartPoint = [1 1 1];
	                tempOpts.Lower = [0 -Inf 0];
					[ tempFitRes, tempGof ] = computeEquationFit( tempQstepAllQP, tempRateAllQP, equationRLib, tempOpts, tempFlagShow );
					tempFitResAllLibPic(idxLibPic, 5:8) = [tempGof.rsquare tempFitRes.u tempFitRes.b tempFitRes.m];
	            end

				% collect data.
				tempFitResLibAllSeq{idxSeq} = tempFitResAllLibPic;
	        end
	        fitResLibAllSeq{idxReso} = tempFitResLibAllSeq;
	    end

        if flagReadSeqListTrain
            fitResLibAllSeq_Train = fitResLibAllSeq;
        elseif flagReadSeqListTest
            fitResLibAllSeq_Test = fitResLibAllSeq;
        end
	end

	%% reference models are verified here.
	if subFlagVerifyReferenceModel
		% prepare data.
		% get basic data.
		if flagReadSeqListTrain
            cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Train;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Train;
        elseif flagReadSeqListTest
            cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
        end
        % collect data.
        fitResKeyAllSeq = cell(numReso, 1);
        for idxReso = 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempFitResKeyAllSeq = cell(tempNumSeq, 1);
			for idxSeq = 1: tempNumSeq
				% prepare data.
	            tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempListLibFreq = cellLibInfo{idxReso}{idxSeq, 3};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
                tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
	            tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length',tempEstInfoAllLibPic)>0);
				tempEstInfoAllKeyPic = oriContentInfoAllSeq{idxReso}{idxSeq, 1};

				% collect data. There are 3 models to test.
				tempFitResAllKeyPic = []; % one for rate and one for distortion.

				for idxLibPic = 1: tempNumLibPic
					tempCountKeyPicPerLibPic = 0;
					for idxKeyPic = 1: tempNumKeyPic
						% prepare data.
						if tempListRefLib(idxKeyPic) + 1 ~= idxLibPic
							continue;
						end
						% skip the key picture which reference itself.
						if tempListOrgLib(idxLibPic) == tempListKeyPic(idxKeyPic)
							tempCountKeyPicPerLibPic = tempCountKeyPicPerLibPic + 1;
	                    	continue;
                        end
						% prepare estimated content information.
				        tempEstResiEnergyLib = tempEstInfoAllLibPic{idxLibPic, 1}{1};
				        tempEstResiEnergyKeyInter = tempEstInfoAllKeyPic{idxKeyPic, 1}{1};
				        tempEstResiEnergyKeyIntra = tempEstInfoAllKeyPic{idxKeyPic, 2}{1};
				        tempPercentInterKey = tempEstInfoAllKeyPic{idxKeyPic, 3};
						% collect data.
						tempFitResEachKeyPic = zeros(1, 10);
						tempCountKeyPicPerLibPic = tempCountKeyPicPerLibPic + 1;

						% loop for each qp.
						for idxQP = 1: numQP
							% collect data.
							tempDistLib = zeros(numDqp, 1);
							tempRateLib = zeros(numDqp, 1);
							tempDistKey = zeros(numDqp, 1);
							tempRateKey = zeros(numDqp, 1);
							for idxDqp = 1: numDqp
								tempRateLib(idxDqp) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(4);
								tempDistLib(idxDqp) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(5);
								tempRateKey(idxDqp) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(tempCountKeyPicPerLibPic, 4);
								tempDistKey(idxDqp) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(tempCountKeyPicPerLibPic, 5);
							end
							% convert data.
							tempDistLib = psnr2mse(tempDistLib);
							tempRateLib = tempRateLib / tempWidth / tempHeight;
							tempDistKey = psnr2mse(tempDistKey);
							tempRateKey = tempRateKey / tempWidth / tempHeight;

                            tempFitResEachKeyPic(1:4) = [tempEstResiEnergyLib tempEstResiEnergyKeyInter tempEstResiEnergyKeyIntra tempPercentInterKey];
							% fit distortion of key pictures with difference reference Q.
							tempFlagShow = 0;
							if tempFlagShow
			                	subplot(2,1,1);
			                end
			                eqDKeyDRef = @(a,b,x) a.*x+b;
		        			equationDKey = eqDKeyDRef;
		                    tempOpts = fitoptions(equationDKey);
		                    tempOpts.StartPoint = [1 1];
		                    tempOpts.Lower = [0 -Inf];
							[ tempFitRes, tempGof ] = computeEquationFit( tempDistLib, tempDistKey, equationDKey, tempOpts, tempFlagShow );
							tempFitResEachKeyPic(5:7) = [tempGof.rsquare tempFitRes.a tempFitRes.b];
	                        if tempFlagShow
	                            hold on;
	                        end

							% fit rate of key pictures with difference reference Q.
							if tempFlagShow
			                	subplot(2,1,2);
			                end
	        				eqRKeyDRef = @(a,b,x) a.*x+b;
		        			equationRKey = eqRKeyDRef;
		                    tempOpts = fitoptions(equationRKey);
		                    tempOpts.StartPoint = [1 1];
		                    tempOpts.Lower = [0 -Inf];
							[ tempFitRes, tempGof ] = computeEquationFit( tempDistLib, tempRateKey, equationRKey, tempOpts, tempFlagShow );
							tempFitResEachKeyPic(8:10) = [tempGof.rsquare tempFitRes.a tempFitRes.b];
	                        if tempFlagShow
	                            hold on;
	                        end

							% collect data.
							tempFitResAllKeyPic = [tempFitResAllKeyPic; tempFitResEachKeyPic];
	                    end
	                    close all;
					end
				end
				% collect data.
				tempFitResKeyAllSeq(idxSeq) = {tempFitResAllKeyPic};
	        end
	        % collect data.
	        fitResKeyAllSeq(idxReso) = {tempFitResKeyAllSeq};
	    end

        if flagReadSeqListTrain
            fitResKeyAllSeq_Train = fitResKeyAllSeq;
        elseif flagReadSeqListTest
            fitResKeyAllSeq_Test = fitResKeyAllSeq;
        end
	end

	%% all model are verified here.
	if subFlagPlotActualRateDist
		% switches.
		subsubFlagCollectData = 1;

        % prepare data.
		% get basic data.
		if flagReadSeqListTrain
            cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Train;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Train;
        elseif flagReadSeqListTest
            cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
        end

        % collect data.
        if subsubFlagCollectData
	        DLibQ = [];
	        RLibQ = [];
	        DLibResi = [];
	        RLibResi = [];
	        DKeyRefQ = [];
	        RKeyRefQ = [];
	        DKeyResiInter = [];
	        RKeyResiInter = [];
	        DKeyResiSim = [];
	        RKeyResiSim = [];

	        for idxReso = 1: numReso
	        	if numAllSeq(idxReso) < 1
					continue;
				end
				% prepare data.
				tempNumSeq = numAllSeq(idxReso);
				for idxSeq = 1: tempNumSeq
					% prepare data.
		            tempWidth = cellListSeq{idxReso}{idxSeq, 5};
		            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
					tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
					tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
					tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
					tempNumLibPic = length(tempListOrgLib);
					tempNumKeyPic = length(tempListKeyPic);
					tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
		            tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length',tempEstInfoAllLibPic)>0);
					tempEstInfoAllKeyPic = oriContentInfoAllSeq{idxReso}{idxSeq, 1};

					% plot curve for each library picture.
					for idxLibPic = 1: tempNumLibPic
						% collect library data.
						tempR = tempEstInfoAllLibPic{idxLibPic}{1};
						tempD = tempEstInfoAllLibPic{idxLibPic}{1};
						for idxQP = 1: numQP
							for idxDqp = 1: numDqp
								tempR = [tempR qp2qstep(listQP(idxQP)+listDqp(idxDqp)) bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(4)/tempWidth/tempHeight];
								tempD = [tempD qp2qstep(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(3)) psnr2mse(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(5))];
							end
						end
						DLibQ = [DLibQ; tempD];
						RLibQ = [RLibQ; tempR];

						% collect key picture data.
						% prepare data.
						tempNumKeyPicPerLib = size(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{1}{1}{idxLibPic, 2}, 1);
						tempList = (tempListRefLib == (idxLibPic - 1));
	                    tempList = find(tempList == 1, tempNumKeyPicPerLib);
						tempEstInfoAllKeyPicPerLibPic = tempEstInfoAllKeyPic(tempList, :);
						% for each key picture.
						for idxKeyPic = 1: tempNumKeyPicPerLib
	                        if ismember(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{1}{1}{idxLibPic, 2}(idxKeyPic, 1), tempListOrgLib)
	                            continue;
	                        end

							tempFlagShow = 0;
							if tempFlagShow
								tempDAllQP = zeros(numQP, numDqp, 2);
								tempRAllQP = zeros(numQP, numDqp, 2);
							end

							for idxQP = 1: numQP
								tempR = [tempEstInfoAllKeyPicPerLibPic{idxKeyPic, 1}{1} tempEstInfoAllKeyPicPerLibPic{idxKeyPic, 2}{1} estInfoAllLibPic{idxLibPic}{1} tempEstInfoAllKeyPicPerLibPic{idxKeyPic, 3} qp2qstep(listQP(idxQP))];
								tempD = [tempEstInfoAllKeyPicPerLibPic{idxKeyPic, 1}{1} tempEstInfoAllKeyPicPerLibPic{idxKeyPic, 2}{1} estInfoAllLibPic{idxLibPic}{1} tempEstInfoAllKeyPicPerLibPic{idxKeyPic, 3} qp2qstep(listQP(idxQP))];
								for idxDqp = 1: numDqp
									tempR = [tempR qp2qstep(listQP(idxQP)+listDqp(idxDqp)) bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(idxKeyPic, 4)/tempWidth/tempHeight];
									tempD = [tempD qp2qstep(listQP(idxQP)+listDqp(idxDqp)) psnr2mse(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(idxKeyPic, 5))];
									if tempFlagShow
										tempDAllQP(idxQP, idxDqp, 1) = psnr2mse(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(idxKeyPic, 5));
										tempRAllQP(idxQP, idxDqp, 1) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(idxKeyPic, 4)/tempWidth/tempHeight;
	                                    tempDAllQP(idxQP, idxDqp, 2) = psnr2mse(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(5));
	                                    tempRAllQP(idxQP, idxDqp, 2) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(4)/tempWidth/tempHeight;
									end
								end
								DKeyRefQ = [DKeyRefQ; tempD];
								RKeyRefQ = [RKeyRefQ; tempR];
							end

	                        if tempFlagShow
	                        	subplot(2,1,1);
	                        	for idxQP = 1: numQP
	                        		plot(tempDAllQP(idxQP, :, 2), tempDAllQP(idxQP, :, 1), '-d');
	                        		hold on;
	                        	end
	                        	title(['distortion of ' cellListSeq{idxReso}{idxSeq, 2}]);
	                        	xlabel('reference distortion');
	                        	ylabel('distortion(MSE)');
	                        	hold off;
	                        	subplot(2,1,2);
	                        	for idxQP = 1: numQP
	                        		plot(tempRAllQP(idxQP, :, 2), tempRAllQP(idxQP, :, 1), '-d');
	                        		hold on;
	                        	end
	                        	title(['rate of ' cellListSeq{idxReso}{idxSeq, 2}]);
	                        	xlabel('reference rate');
	                        	ylabel('rate(bpp)');
	                        	hold off;
	                        end
	                    end
					end
	            end
	        end
            if flagReadSeqListTrain
                DLibQ_Train = DLibQ;
                RLibQ_Train = RLibQ;
                DLibResi_Train = DLibResi;
                RLibResi_Train = RLibResi;
                DKeyRefQ_Train = DKeyRefQ;
                RKeyRefQ_Train = RKeyRefQ;
                DKeyResiInter_Train = DKeyResiInter;
                RKeyResiInter_Train = RKeyResiInter;
                DKeyResiSim_Train = DKeyResiSim;
                RKeyResiSim_Train = RKeyResiSim;
            elseif flagReadSeqListTest
                DLibQ_Test = DLibQ;
                RLibQ_Test = RLibQ;
                DLibResi_Test = DLibResi;
                RLibResi_Test = RLibResi;
                DKeyRefQ_Test = DKeyRefQ;
                RKeyRefQ_Test = RKeyRefQ;
                DKeyResiInter_Test = DKeyResiInter;
                RKeyResiInter_Test = RKeyResiInter;
                DKeyResiSim_Test = DKeyResiSim;
                RKeyResiSim_Test = RKeyResiSim;
            end
            

			% prepare data for library rate/distortion-resi.
			tempResi = DLibQ(:, 1);
			tempQ = DLibQ(1, 2:2:end);
			tempNumCol = size(DLibQ, 1) * 2 + 1;
			tempNumRow = length(tempQ);
			DLibResi = zeros(tempNumRow, tempNumCol);
			RLibResi = zeros(tempNumRow, tempNumCol);
			DLibResi(:, 1) = tempQ;
			RLibResi(:, 1) = tempQ;
			DLibResi(:, 2:2:end) = ones(tempNumRow, 1) * tempResi';
			RLibResi(:, 2:2:end) = ones(tempNumRow, 1) * tempResi';
			DLibResi(:, 3:2:end) = DLibQ(:, 3:2:end)';
			RLibResi(:, 3:2:end) = RLibQ(:, 3:2:end)';

			% prepare data for key rate/distortion-resi.
			for idxQP = 1: numQP
				tempResi = DKeyRefQ(idxQP:numQP:end, 1);
				tempRefQ = DKeyRefQ(idxQP, 6:2:end);
				tempNumCol = size(DKeyRefQ, 1) / numQP * 2 + 1;
				tempNumRow = length(tempRefQ);
				tempDKeyRefQ = zeros(tempNumRow, tempNumCol);
				tempRKeyRefQ = zeros(tempNumRow, tempNumCol);
				tempDKeyRefQ(:, 1) = tempRefQ;
				tempRKeyRefQ(:, 1) = tempRefQ;
				tempDKeyRefQ(:, 2:2:end) = ones(tempNumRow, 1) * tempResi';
				tempRKeyRefQ(:, 2:2:end) = ones(tempNumRow, 1) * tempResi';
				tempDKeyRefQ(:, 3:2:end) = DKeyRefQ(idxQP:numQP:end, 7:2:end)';
				tempRKeyRefQ(:, 3:2:end) = RKeyRefQ(idxQP:numQP:end, 7:2:end)';
				DKeyResiInter = [DKeyResiInter; tempDKeyRefQ];
				RKeyResiInter = [RKeyResiInter; tempRKeyRefQ];
			end

			% prepare data for key rate/distortion-sim.
			for idxQP = 1: numQP
				tempSim = DKeyRefQ(idxQP:numQP:end, 3);
				tempRefQ = DKeyRefQ(idxQP, 6:2:end);
				tempNumCol = size(DKeyRefQ, 1) / numQP * 2 + 1;
				tempNumRow = length(tempRefQ);
				tempDKeyRefQ = zeros(tempNumRow, tempNumCol);
				tempRKeyRefQ = zeros(tempNumRow, tempNumCol);
				tempDKeyRefQ(:, 1) = tempRefQ;
				tempRKeyRefQ(:, 1) = tempRefQ;
				tempDKeyRefQ(:, 2:2:end) = ones(tempNumRow, 1) * tempSim';
				tempRKeyRefQ(:, 2:2:end) = ones(tempNumRow, 1) * tempSim';
				tempDKeyRefQ(:, 3:2:end) = DKeyRefQ(idxQP:numQP:end, 7:2:end)';
				tempRKeyRefQ(:, 3:2:end) = RKeyRefQ(idxQP:numQP:end, 7:2:end)';
				DKeyResiSim = [DKeyResiSim; tempDKeyRefQ];
				RKeyResiSim = [RKeyResiSim; tempRKeyRefQ];
			end
		end

		% plot data.
        % common model.
        eqDLib = @(a,e,q) a.*e.^0.5.*q;
		eqRLib = @(u,f,b,e,q) u.*exp(-f.*q./(e.^0.5))+b;
		eqDKey = @(a1,a2,e1,e2,e3,s,q1,q2) (1-s).*eqDLib(a2,e1,q1)+s.*eqDLib(a1,e2+eqDLib(a2,e3,q2),q1);
		eqRKey = @(u1,f1,b1,u2,f2,b2,a2,e1,e2,e3,s,q1,q2) (1-s).*eqRLib(u2,f2,b2,e1,q1)+s.*eqRLib(u1,f1,b1,e2+eqDLib(a2,e3,q2),q1);
		% constant coefficients.
		a1=0.1216;a2=0.08302;
		u1=3.578;f1=0.947;b1=0;
		u2=4.645;f2=1.906;b2=0.2732;
		wD=1.647;wR=0.155;
		% plot library data.
		subplot(2,2,1);
		tempNumRow = size(DLibQ, 1);
		for idx = 1: tempNumRow
			plot(DLibQ(idx, 2:2:end), DLibQ(idx, 3:2:end), 'x');
            hold on;
            plot(DLibQ(idx, 2:2:end), eqDLib(a2,DLibQ(idx, 1),DLibQ(idx, 2:2:end)), 'rx');
			hold off;
		end
		hold off;
		subplot(2,2,2);
		tempNumRow = size(RLibQ, 1);
		for idx = 1: tempNumRow
			plot(RLibQ(idx, 2:2:end), RLibQ(idx, 3:2:end), 'x');
			hold on;
            plot(RLibQ(idx, 2:2:end), eqRLib(u2,f2,b2,RLibQ(idx, 1),RLibQ(idx, 2:2:end)), 'rx');
			hold off;
		end
		hold off;
		subplot(2,2,3);
		tempNumRow = size(DLibResi, 1);
		for idx = 1: tempNumRow
			plot(DLibResi(idx, 2:2:end), DLibResi(idx, 3:2:end), 'x');
			hold on;
		end
		hold off;
		subplot(2,2,4);
		tempNumRow = size(RLibResi, 1);
		for idx = 1: tempNumRow
			plot(RLibResi(idx, 2:2:end), RLibResi(idx, 3:2:end), 'x');
			hold on;
		end
		hold off;
		% plot key data.
        subplot(3,2,1);
		tempNumRow = size(DKeyRefQ, 1);
		for idx = 1: tempNumRow
			plot(DKeyRefQ(idx, 6:2:end), DKeyRefQ(idx, 7:2:end), 'x');
            hold on;
            e1=DKeyRefQ(idx, 2);e2=DKeyRefQ(idx, 1);e3=DKeyRefQ(idx, 3);
            s=DKeyRefQ(idx, 4);q1=DKeyRefQ(idx, 5);
            plot(DKeyRefQ(idx, 6:2:end), eqDKey(a1,a2,e1,e2,e3,s,q1,DKeyRefQ(idx, 6:2:end)), 'rx');
			hold off;
		end
		hold off;
		subplot(3,2,2);
		tempNumRow = size(RKeyRefQ, 1);
		for idx = 1: tempNumRow
			plot(RKeyRefQ(idx, 6:2:end), RKeyRefQ(idx, 7:2:end), 'x');
			hold on;
            e1=RKeyRefQ(idx, 2);e2=RKeyRefQ(idx, 1);e3=RKeyRefQ(idx, 3);
            s=RKeyRefQ(idx, 4);q1=RKeyRefQ(idx, 5);
            plot(RKeyRefQ(idx, 6:2:end), eqRKey(u1,f1,b1,u2,f2,b2,a2,e1,e2,e3,s,q1,RKeyRefQ(idx, 6:2:end)), 'rx');
			hold off;
		end
		hold off;
		subplot(3,2,3);
		tempNumRow = size(DKeyResiInter, 1);
		for idx = 1: tempNumRow
			plot(DKeyResiInter(idx, 2:2:end), DKeyResiInter(idx, 3:2:end), 'x');
			hold on;
		end
		hold off;
		subplot(3,2,4);
		tempNumRow = size(RKeyResiInter, 1);
		for idx = 1: tempNumRow
			plot(RKeyResiInter(idx, 2:2:end), RKeyResiInter(idx, 3:2:end), 'x');
			hold on;
		end
		hold off;
		subplot(3,2,5);
		tempNumRow = size(DKeyResiSim, 1);
		for idx = 1: tempNumRow
			plot(DKeyResiSim(idx, 2:2:end), DKeyResiSim(idx, 3:2:end), 'x');
			hold on;
		end
		hold off;
		subplot(3,2,6);
		tempNumRow = size(RKeyResiSim, 1);
		for idx = 1: tempNumRow
			plot(RKeyResiSim(idx, 2:2:end), RKeyResiSim(idx, 3:2:end), 'x');
			hold on;
		end
		hold off;
        
	end

	%% plot the trend curve of rate and distortion based on theory.
	if subFlagPlotTheoryRateDistTrend
        % fitting model.
        equationIntegral = @(a,b,x) (-0.5.*exp(-a.*x).*((a.*x+1).^2+1)+1+1./(1-exp(-x)).*(-0.5*exp(-(1+a).*x)).*(((1-b).*x+1).^2+1-exp(x).*((b.*x-1).^2+1)));
        equationDLib = @(a,b,c,x) c.*equationIntegral(a,b,x./((c./2).^0.5));
        equationRLib = @(a,b,c,d,e,x) (-d*log(equationIntegral(a,b,x./((c./2).^0.5)))+e);
        equationEntropy = @(a,x) (exp(-a.*x).*(log2((1-exp(-a.*x))./(1-exp(-x)))+1+x.*(a.*(1-exp(-x))+exp(-x))./(1-exp(-x))/log(2))-log2(1-exp(-a.*x)));
        equationREntropy = @(a,c,f,g,x) (f*equationEntropy(a,x./((c./2).^0.5))+g);
        aIntra = 2/3; bIntra = 1-aIntra;
        aInter = 5/6; bInter = 1-aInter;

        % prepare data.
        qList = qp2qstep([16:2:40]');
        resiEner = [ 20 50 150 250 400 500]';

        % show rate distortion with Q self.
        equationR = @(c,x) (equationEntropy(aIntra,x./((c./2).^0.5)));
		equationD = @(c,x) c.*equationIntegral(aIntra,bIntra,x./((c./2).^0.5));
        equationRSimp = @(c,x) exp(-x./((c./2).^0.5));
        equationDSimp = @(c,x) c.^0.5.*x;
        % plot rate trend with q self.
        subplot(3,2,1);
        data = [];
        for idx = 1: length(resiEner)
        	plot(qList, log(equationR(resiEner(idx), qList)), '-d');
            %plot(qList, log(equationRSimp(resiEner(idx), qList)), '-rd');
        	hold on;
        end
        title('rate under different residue energy');
        xlabel('Q'); ylabel('log rate(bpp)');
        hold off;
        % plot distortion with q self.
        subplot(3,2,2);
        for idx = 1: length(resiEner)
        	plot(qList, equationD(resiEner(idx), qList), '-d');
            %plot(qList, equationDSimp(resiEner(idx), qList), '-rd');
            %data = [data equationD(resiEner(idx), qList)];
        	hold on;
        end
        title('distortion under different residue energy');
        xlabel('Q'); ylabel('distortion(MSE per pixel)');
        hold off;
        % plot rate trend with residue energy.
        subplot(3,2,3);
        for idx = 1: length(qList)
        	plot(resiEner.^(-0.5), log(equationR(resiEner, qList(idx))), '-d');
            %plot(log(resiEner), equationRSimp(resiEner, qList(idx)), '-rd');
        	hold on;
        end
        title('rate under difference Q');
        xlabel('resiEnergy.(-0.5)'); ylabel('log rate(bpp)');
        hold off;
        % plot distortion with residue energy.
        subplot(3,2,4);
        for idx = 1: length(qList)
        	plot(resiEner.^(0.5), equationD(resiEner, qList(idx)), '-d');
            data = [data equationD(resiEner, qList(idx))];
            %plot(resiEner, equationDSimp(resiEner, qList(idx)), '-rd');
        	hold on;
        end
        title('distortion under different Q');
        xlabel('resiEnergy.(0.5)'); ylabel('distortion(MSE per pixel)');
        hold off;

        % show rate and distortion with self Q and reference Q.
        equationRRef = @(c,d,q,x) (equationEntropy(aInter,x./(((d+c.*equationD(c,q))./2).^0.5)));
        equationDRef = @(c,d,q,x) ((d+c.*equationD(c,q)).*equationIntegral(aInter,bInter,x./(((d+c.*equationD(c,q))./2).^0.5)));
        resiRef = 100;
        diff = 0.1;
        qSelfList = qp2qstep(22:5:37);
        qRefList = qp2qstep(11:1:40);
        % plot rate.
        subplot(3,2,5);
        for idx = 1: length(qSelfList)
        	plot(qRefList, equationRRef(resiRef, diff, qRefList, qSelfList(idx)), '-d');
        	hold on;
        end
        title(['rate under difference self Q with diff=' num2str(diff)]);
        xlabel('reference Q'); ylabel('rate(per pixel)');
        hold off;
        % plot distortion.
        subplot(3,2,6);
        for idx = 1: length(qSelfList)
        	plot(qRefList, equationDRef(resiRef, diff, qRefList, qSelfList(idx)), '-d');
        	hold on;
        end
        title(['distortion under difference self Q with diff=' num2str(diff)]);
        xlabel('reference Q'); ylabel('distortion(MSE per pixel)');
        hold off;
	end

	%% train the model coefficients.
	if subFlagTrainModelCoeffLib
		% prepare data.
		if flagReadSeqListTrain
            cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Train;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Train;
        elseif flagReadSeqListTest
            cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
            cellEncAnalInfo = cellEncAnalInfo_Test;
        end
        % switches.
        subsubFlagCollectData = 1;
		% collect data.
		fitResLibAllSeq = cell(numReso, 1);
		trainDataLib = cell(numReso, 2); % for 4 train experiments.
        for idxReso = 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempFitResLibAllSeq = cell(tempNumSeq, 1);
			tempTrainDataDLib = [];
			tempTrainDataRLib = [];
			for idxSeq = 1: tempNumSeq
				% prepare data.
	            tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempListLibFreq = cellLibInfo{idxReso}{idxSeq, 3};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
				tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
	            tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length',tempEstInfoAllLibPic)>0);
				tempEstInfoAllKeyPic = oriContentInfoAllSeq{idxReso}{idxSeq, 1};

	            % collect data. There are 3 models to test.
				tempFitResAllLibPic = zeros(tempNumLibPic, 10); % for distortion and rate.
				if subsubFlagCollectData
					% collect library picture data.
					for idxLibPic = 1: tempNumLibPic
						% collect data.
						tempDistAllQP = zeros(numQP * numDqp, 1);
						tempRateAllQP = zeros(numQP * numDqp, 1);
						tempQstepAllQP = zeros(numQP * numDqp, 1);

						% prepare estimated content information.
			        	%tempEstResiEnergyLib = tempEstInfoAllLibPic{idxLibPic, 1}{1};
                        tempEstResiEnergyLib = cellEncAnalInfo{idxReso}{idxSeq}{4}{1}{1}(idxLibPic, 3);
						for idxQP = 1: numQP
							for idxDqp = 1: numDqp
								tempQstepAllQP((idxQP-1)*numDqp + idxDqp, 1) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(3);
								tempRateAllQP((idxQP-1)*numDqp + idxDqp, 1) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(4);
								tempDistAllQP((idxQP-1)*numDqp + idxDqp, 1) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(5);
							end
						end
						% remove repeat data.
						[tempQstepAllQP, idxUnique, ~] = unique(tempQstepAllQP, 'rows');
		                tempQstepAllQP = qp2qstep(tempQstepAllQP);
						tempRateAllQP = tempRateAllQP(idxUnique) / tempWidth / tempHeight;
						tempDistAllQP = psnr2mse(tempDistAllQP(idxUnique));
						% fit data.
						tempFlagShow = 1;
	                    if tempFlagShow
	                        subplot(2,1,1);
	                    end
						% fit rate model of library pictures.
	        			%eqDLib = @(a,x) a.*tempEstResiEnergyLib.^0.5.*x;
                        eqDLib = @(a,x) a.*tempEstResiEnergyLib.*x;
		                tempOpts = fitoptions(eqDLib);
		                tempOpts.StartPoint = [1];
		                tempOpts.Lower = [1e-10];
						[ tempFitRes, tempGof ] = computeEquationFit( tempQstepAllQP, tempDistAllQP, eqDLib, tempOpts, tempFlagShow );
						tempFitResAllLibPic(idxLibPic, 1:5) = [tempGof.rsquare tempFitRes.a 0 0 0];
                        %if tempGof.rsquare > 0.9
                            tempTrainDataDLib = [tempTrainDataDLib; tempEstResiEnergyLib*ones(size(tempQstepAllQP)) tempQstepAllQP tempDistAllQP];
                        %end
                        if tempFlagShow
                            hold on;
                        end

	                    if tempFlagShow
	                        subplot(2,1,2);
	                    end
		        		eqRLib = @(u,b,x) u.*tempEstResiEnergyLib.^0.5./x+b;
		                tempOpts = fitoptions(eqRLib);
		                tempOpts.StartPoint = [1 1];
		                tempOpts.Lower = [0 -Inf];
						[ tempFitRes, tempGof ] = computeEquationFit( tempQstepAllQP, tempRateAllQP, eqRLib, tempOpts, tempFlagShow );
						tempFitResAllLibPic(idxLibPic, 6:10) = [tempGof.rsquare tempFitRes.u tempFitRes.b 0 0];
                        tempRateAllQP = tempRateAllQP - tempFitRes.b;
                        %if tempGof.rsquare > 0.9
                            tempTrainDataRLib = [tempTrainDataRLib; tempEstResiEnergyLib*ones(size(tempQstepAllQP)) tempQstepAllQP tempRateAllQP];
                        %end
                        if tempFlagShow
                            hold on;
                        end
					end
	            end
	            
				% collect data.
				tempFitResLibAllSeq(idxSeq) = {tempFitResAllLibPic};
	        end
	        % collect data.
	        fitResLibAllSeq(idxReso) = {tempFitResLibAllSeq};
			trainDataLib(idxReso, :) = [{tempTrainDataDLib} {tempTrainDataRLib}];
	    end

        if flagReadSeqListTrain && subsubFlagCollectData
            fitResLibAllSeq_Train = fitResLibAllSeq;
            trainDataLib_Train = trainDataLib;
        elseif flagReadSeqListTest && subsubFlagCollectData
            fitResLibAllSeq_Test = fitResLibAllSeq;
            trainDataLib_Test = trainDataLib;
        end
        clear temp* idx*;

        % train general model for all sequences.
        if flagReadSeqListTrain
            trainDataLib = trainDataLib_Train;
        elseif flagReadSeqListTest
            trainDataLib = trainDataLib_Test;
        end
    	for idxReso = 1: numReso
    		if numAllSeq(idxReso) < 1
				continue;
            end
            tempTrainData = trainDataLib{idxReso, 1};
            % train model for distortion.
            tempResiEnerList = tempTrainData(:, 1);
            tempQList = tempTrainData(:, 2);
            tempDistList = tempTrainData(:, 3);
            %eqDLibSimp = @(a,c,x) a.*c.^0.5.*x;
            eqDLibSimp = @(a,c,x) a.*c.*x;
            tempFt = fittype(eqDLibSimp, 'independent', {'c','x'}, 'dependent', 'z', 'coefficients', 'a');
            tempOpts = fitoptions(eqDLibSimp);
            tempOpts.StartPoint = [1];
            tempOpts.Lower = [0];
            [tempFitRes, tempGof] = fit([tempResiEnerList tempQList], tempDistList, tempFt, tempOpts);
            tempFitRes
            tempGof
            plot( tempFitRes, [tempResiEnerList tempQList], tempDistList );

            % train model for rate.
            tempTrainData = trainDataLib{idxReso, 2};
            tempResiEnerList = tempTrainData(:, 1);
            tempQList = tempTrainData(:, 2);
            tempRateList = tempTrainData(:, 3);
            equationRLibSimp = @(u,c,x) u.*c.^0.5./x;
            tempFt = fittype(equationRLibSimp, 'independent', {'c','x'}, 'dependent', 'z', 'coefficients', {'u'});
            tempOpts = fitoptions(equationRLibSimp);
            tempOpts.StartPoint = [1];
            tempOpts.Lower = [0];
            [tempFitRes, tempGof] = fit([tempResiEnerList tempQList], tempRateList, tempFt, tempOpts);
            tempFitRes
            tempGof
            plot( tempFitRes, [tempResiEnerList tempQList], tempRateList );
        end
	end

	%% train the model coefficients.
	if subFlagTrainModelCoeffKey
		% prepare data.
		if flagReadSeqListTrain
            cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Train;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Train;
            cellEncAnalInfo = cellEncAnalInfo_Train;
        elseif flagReadSeqListTest
            cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
            bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
            cellEncAnalInfo = cellEncAnalInfo_Test;
        end
        % switches.
        subsubFlagCollectData = 1;
		% collect data.
		fitResKeyAllSeq = cell(numReso, 1); % one for library picture and one for key pictures.
		trainDataKey = cell(numReso, 2); % for 4 train experiments.
        sumScale = [];
        for idxReso = 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
			% collect data.
			tempFitResKeyAllSeq = cell(tempNumSeq, 2);
			tempTrainDataKeyDKey = [];
			tempTrainDataKeyRKey = [];
			for idxSeq = 1: tempNumSeq
				% prepare data.
	            tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempListLibFreq = cellLibInfo{idxReso}{idxSeq, 3};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
				tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
	            tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length',tempEstInfoAllLibPic)>0);
				tempEstInfoAllKeyPic = oriContentInfoAllSeq{idxReso}{idxSeq, 1};

	            % collect data.
				tempFitResAllKeyPic = zeros(tempNumKeyPic*numQP, 10); % one for rate and one for distortion.
				tempFitResScaleDistEAllKeyPic = zeros(tempNumKeyPic*numQP, 10); % one for rate and one for distortion.
				tempFitResRateAllKeyPic = zeros(tempNumKeyPic*numQP, 10);
                if subsubFlagCollectData 
					% collect key picture data.
					tempCountKeyPicPerLib = zeros(tempNumLibPic, 1) + 1;
                    tempSumScale = zeros(tempNumLibPic, numQP);
					for idxKeyPic = 1: tempNumKeyPic
						% prepare data.
						idxLibPic = tempListRefLib(idxKeyPic) + 1;
						idxKeyPicPerLib = tempCountKeyPicPerLib(idxLibPic);
						% skip the key picture which reference itself.
						if tempListOrgLib(tempListRefLib(idxKeyPic) + 1) == tempListKeyPic(idxKeyPic)
							tempCountKeyPicPerLib(idxLibPic) = tempCountKeyPicPerLib(idxLibPic) + 1;
		                    continue;
                        end
						% prepare estimated content information.
				        tempEstResiEnergyLib = tempEstInfoAllLibPic{idxLibPic, 1}{1};
				        tempEstResiEnergyKeyInter = tempEstInfoAllKeyPic{idxKeyPic, 1}{1};
				        tempEstResiEnergyKeyIntra = tempEstInfoAllKeyPic{idxKeyPic, 2}{1};
				        tempPercentInterKey = tempEstInfoAllKeyPic{idxKeyPic, 3};
                        
						% loop for each qp.
						for idxQP = 1: numQP
		                    % prepare data.
		                    tempQpKey = listQP(idxQP);
		                    tempQKey = qp2qstep(tempQpKey);
                            tempEncInterCount = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 4);
                            tempPercentInterKey = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 2);
                            tempPercentSKipKey = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 3);
                            if tempEncInterCount ~= 0
                                tempEstResiEnergyKeyInter = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 5) * tempWidth * tempHeight / tempEncInterCount;
                            else
                                tempEstResiEnergyKeyInter = 0;
                            end
                            
							% collect data.
							tempDistLib = zeros(numDqp, 1);
							tempDistKey = zeros(numDqp, 1);
							tempRateKey = zeros(numDqp, 1);
                            tempPercent = zeros(numDqp, 2);
                            tempEstEKeyInter = zeros(numDqp, 1);
                            tempEstEKeyIntra = zeros(numDqp, 1);
							for idxDqp = 1: numDqp
								tempDistLib(idxDqp) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(5);
								tempRateKey(idxDqp) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(idxKeyPicPerLib, 4);
								tempDistKey(idxDqp) = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(idxKeyPicPerLib, 5);
                                tempPercent(idxDqp, 1) = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 2);
                                tempPercent(idxDqp, 2) = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 3);
                                tempEstEKeyInter(idxDqp) = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 5);
								tempEstEKeyIntra(idxDqp) = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 8);
                            end
                            tempPercentInterKeyMean = mean(tempPercent(:, 1));
                            tempPercentSkipKeyMean = mean(tempPercent(:, 2));
							% convert data.
							tempDistLib = psnr2mse(tempDistLib);
							tempDistKey = psnr2mse(tempDistKey);
							tempRateKey = tempRateKey / tempWidth / tempHeight;

							% fit distortion of key pictures with reference distortion.
							tempFlagShow = 1;
							if tempFlagShow
			                	subplot(2,1,1);
			                end
			                eqDKeyDRef = @(a,b,x) a.*x+b;
		        			equationDKey = eqDKeyDRef;
		                    tempOpts = fitoptions(equationDKey);
		                    tempOpts.StartPoint = [1 1];
		                    tempOpts.Lower = [0 -Inf];
							[ tempFitRes, tempGof ] = computeEquationFit( tempDistLib, tempDistKey, equationDKey, tempOpts, tempFlagShow );
							tempFitResAllKeyPic((idxKeyPic-1)*numQP+idxQP, 1:10) = [tempGof.rsquare tempFitRes.a tempEstResiEnergyKeyInter tempPercentInterKey tempPercentSKipKey tempPercentInterKeyMean tempPercentSkipKeyMean tempQKey 0 0];
	                        if tempFlagShow
	                            hold on;
                            end
                            if tempGof.rsquare > 0.9
                                tempTrainDataKeyDKey = [tempTrainDataKeyDKey; ones(size(tempDistLib))*[tempEstResiEnergyKeyInter;tempPercentInterKey;tempQKey;]' tempDistLib (tempDistKey - tempFitRes.b)];
                            end
                            tempSumScale(idxLibPic, idxQP) = tempSumScale(idxLibPic, idxQP) + 0.011*tempPercentInterKey*tempQKey;

							% fit rate of key pictures with difference reference Q.
							if tempFlagShow
			                	subplot(2,1,2);
			                end
	        				eqRKeyDRef = @(a,b,x) a.*x+b;
		        			equationRKey = eqRKeyDRef;
		                    tempOpts = fitoptions(equationRKey);
		                    tempOpts.StartPoint = [1 1];
		                    tempOpts.Lower = [0 -Inf];
							[ tempFitRes, tempGof ] = computeEquationFit( tempDistLib, tempRateKey, equationRKey, tempOpts, tempFlagShow );
							tempFitResRateAllKeyPic((idxKeyPic-1)*numQP+idxQP, 1:10) = [tempGof.rsquare tempFitRes.a tempFitRes.b 0 0 0 0 0 0 0];
	                        if tempFlagShow
	                            hold on;
                            end
                            tempRateKey = tempRateKey - tempFitRes.b;
                            if tempGof.rsquare > 0.9
                                tempTrainDataKeyRKey = [tempTrainDataKeyRKey; ones(size(tempDistLib))*[tempEstResiEnergyKeyIntra;tempEstResiEnergyKeyInter;tempEstResiEnergyLib;tempPercentInterKey;tempQKey;]' tempDistLib tempRateKey];
                            end

% 							% fit key distortion with residue energy to find the scale.
% 							if tempFlagShow
% 								subplot(2,1,2);
% 							end
% 
% 							tempEstELibKey = tempPercent(:, 1) .* tempEstEKeyInter + (1 - tempPercent(:, 1)) .* tempEstEKeyIntra;
% 							eqDKeyEKey = @(a,b,x) a.*x+b;
% 							tempOpts = fitoptions(eqDKeyEKey);
% 							tempOpts.StartPoint = [1 1];
% 							tempOpts.Lower = [0 -Inf];
% 							[ tempFitRes, tempGof ] = fit( tempEstELibKey.*tempQKey, tempDistKey, eqDKeyEKey, tempOpts);
% 							tempFitResScaleDistEAllKeyPic((idxKeyPic-1)*numQP+idxQP, 1:10) = [tempGof.rsquare tempFitRes.a tempEstResiEnergyKeyInter tempPercentInterKey tempPercentSKipKey tempPercentInterKeyMean tempPercentSkipKeyMean tempQKey 0 0];

						end
						% count the number of key pictures.
						tempCountKeyPicPerLib(idxLibPic) = tempCountKeyPicPerLib(idxLibPic) + 1;
					end
				end
				% collect data.
				tempFitResKeyAllSeq(idxSeq, 1) = {tempFitResAllKeyPic};
				tempFitResKeyAllSeq(idxSeq, 2) = {tempFitResScaleDistEAllKeyPic};
                sumScale = [sumScale; tempSumScale];
	        end
	        % collect data.
	        fitResKeyAllSeq(idxReso) = {tempFitResKeyAllSeq};
			trainDataKey(idxReso, :) = [{tempTrainDataKeyDKey} {tempTrainDataKeyRKey}];
	    end

        if flagReadSeqListTrain && subsubFlagCollectData
            fitResKeyAllSeq_Train = fitResKeyAllSeq;
            trainDataKey_Train = trainDataKey;
        elseif flagReadSeqListTest && subsubFlagCollectData
            fitResKeyAllSeq_Test = fitResKeyAllSeq;
            trainDataKey_Test = trainDataKey;
        end
        clear temp* idx*;

        % train general model for all sequences.
        if flagReadSeqListTrain
            trainDataKey = trainDataKey_Train;
        elseif flagReadSeqListTest
            trainDataKey = trainDataKey_Test;
        end
        tempAllTrainData = [];
        for idxReso = 1: numReso
            if numAllSeq(idxReso) < 1
                continue;
            end
            tempAllTrainData = [tempAllTrainData; trainDataKey{idxReso, 1}];
            % train model for distortion.
            tempResiEnerListInter = trainDataKey{idxReso, 1}(:, 1);
            tempPercent = trainDataKey{idxReso, 1}(:, 2);
            tempQKeyList = trainDataKey{idxReso, 1}(:, 3);
            tempDistLib = trainDataKey{idxReso, 1}(:, 4);
            tempDistList = trainDataKey{idxReso, 1}(:, 5);
% 			aIntra = 0.08302;
% 			xList = (tempResiEnerListInter+aIntra*tempResiEnerListLib.^0.5.*qList).^0.5.*tempQKeyList.*tempPercent;
% 			yList = (1-tempPercent).*aIntra.*tempResiEnerListIntra.^0.5.*tempQKeyList;
%             equationDKeySimp = @(a,x,y) a.*x+y;
%         	ft = fittype(equationDKeySimp, 'independent', {'x','y'}, 'dependent', 'z', 'coefficients', 'a');
%         	opts = fitoptions(ft);
%         	opts.StartPoint = [1];
%         	opts.Lower = [1e-10];
%         	[fitRes, gof] = fit([xList yList], tempDistList, ft, opts);
%         	fitRes
%             gof

            xDistList = tempPercent.*tempQKeyList.*tempDistLib;
            tempDistList = tempDistList(xDistList>0);
            xDistList = xDistList(xDistList>0);
            eqDKeySimp = @(a,b,x) a.*x+b;
            tempOpts = fitoptions(eqDKeySimp);
            tempOpts.StartPoint = [1 1];
            tempOpts.Lower = [0 -Inf];
            [tempFitRes, tempGof] = fit(xDistList, tempDistList, eqDKeySimp, tempOpts);
            tempFitRes
            tempGof
            plot( tempFitRes, xDistList, tempDistList );

            % train model for rate.
% 	        	tempResiEnerListIntra = trainData{idxReso, 4}(:, 1);
% 	            tempResiEnerListInter = trainData{idxReso, 4}(:, 2);
% 	            tempResiEnerListLib = trainData{idxReso, 4}(:, 3);
% 	            tempPercent = trainData{idxReso, 4}(:, 4);
% 	            tempQKeyList = trainData{idxReso, 4}(:, 5);
% 	            tempDistLib = trainData{idxReso, 4}(:, 6);
% 	            tempRateList = trainData{idxReso, 4}(:, 7);
%           u2 = 4.645; f2 = 1.906; b2 = 0.2732;
% 			a2 = 0.08302;
% 			xList = tempQKeyList./((tempResiEnerListInter+a2*tempResiEnerListLib.^0.5.*qList).^0.5);
% 			yList = (1./tempPercent-1).*equationRLibSimp(u2,f2,b2,tempResiEnerListIntra,tempQKeyList);
% 			zList = tempRateList ./ tempPercent;
%         	equationRKeySimp = @(u,f,b,x,y) y+u.*exp(-f.*x)+b;
%         	ft = fittype(equationRKeySimp, 'independent', {'x','y'}, 'dependent', 'z', 'coefficients', {'u','f','b'});
%         	tempOpts = fitoptions(ft);
%         	tempOpts.StartPoint = [1 1 1];
%         	tempOpts.Lower = [1e-10 1e-10 0];
%         	[tempFitRes, tempGof] = fit([xList yList], zList, ft, tempOpts);
%             tempFitRes
%             tempGof
% 	            xRateList = tempPercent./tempQKeyList./tempResiEnerListInter.^0.5.*tempDistLib;
% 	            eqRKeySimp = @(a,b,x) a.*x+b;
% 	        	tempOpts = fitoptions(eqRKeySimp);
% 	        	tempOpts.StartPoint = [1 1];
% 	        	tempOpts.Lower = [0 -Inf];
% 	        	[tempFitRes, tempGof] = fit(xRateList, tempRateList, eqRKeySimp, tempOpts);
% 	        	tempFitRes
% 	            tempGof
%                 plot( tempFitRes, xRateList, tempRateList );
        end
        xDistList = tempAllTrainData(:,2).*tempAllTrainData(:,3).*tempAllTrainData(:,4);
        tempDistList = tempAllTrainData(xDistList>0, 5);
        xDistList = xDistList(xDistList>0);
        eqDKeySimp = @(a,b,x) a.*x+b;
        tempOpts = fitoptions(eqDKeySimp);
        tempOpts.StartPoint = [1 1];
        tempOpts.Lower = [0 -Inf];
        [tempFitRes, tempGof] = fit(xDistList, tempDistList, eqDKeySimp, tempOpts);
        tempFitRes
        tempGof
        plot( tempFitRes, xDistList, tempDistList );
	end

	%% train the bit and distortion scale between common pictures and key pictures.
	if subFlagTrainScalePicDivKeyRD
        if flagReadSeqListTrain
            cellListSeq = cellListSeq_Train;
			numAllSeq = numSeq_Train;
			cellLibInfo = cellLibInfo_Train;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Train;
			diffPicAllSeq = diffPicAllSeq_Train;
            bitPsnrAllPicAllLibAllSeq = bitPsnrAllPicAllLibAllSeq_Train;
        elseif flagReadSeqListTest
            cellListSeq = cellListSeq_Test;
			numAllSeq = numSeq_Test;
			cellLibInfo = cellLibInfo_Test;
			oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
			diffPicAllSeq = diffPicAllSeq_Test;
            bitPsnrAllPicAllLibAllSeq = bitPsnrAllPicAllLibAllSeq_Test;
        end

		% show scale in one QP with difference delta QP.
		scaleD = cell(numReso, 1);
		listWAllSeq = cell(numReso, 1);
        distAllLayer =cell(numReso, 1);
		for idxReso = 1: numReso
            if numAllSeq(idxReso) < 1
                continue;
            end
            % prepare data.
            tempNumSeq = numAllSeq(idxReso);
			tempScaleDAllQP = cell(numQP, 1);
            tempDistAllLayer = cell(numQP, 1);
            tempWAllSeq = cell(tempNumSeq, 1);
			for idxQP = 1: numQP
                
	            % collect data.
	            tempScaleAllDQP = [];
	            
				for idxSeq = 1: tempNumSeq
					% prepare data.
					tempListLibFreq = cellLibInfo{idxReso}{idxSeq, 3};
					tempListIntraPic = cellLibInfo{idxReso}{idxSeq, 4};
					tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
					tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
					tempNumLibPic = length(tempListOrgLib);
					tempNumIntraPic = length(tempListIntraPic);
					[~, tempListIdxLibInKeyPic] = ismember(tempListOrgLib, tempListIntraPic);
	                
	                for idxLibPic = 1: tempNumLibPic
                        tempNumKeyPicPerLib = size(bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{1}{idxLibPic},1);
                        for idxKeyPic = 1: tempNumKeyPicPerLib
                        	% prepare data.
                        	tempDiffPic = diffPicAllSeq{idxReso}{idxSeq}{idxLibPic}(idxKeyPic);
                            % collect data.
                            tempDist = zeros(numDqp, 1);
                            tempDistAllLayer = zeros(numDqp, 1);
                            tempList = tempListIntraPic(tempListRefLib == (idxLibPic - 1));
                            tempIdx = find(tempListIntraPic == tempList(idxKeyPic), 1, 'first');
                            if flagReadSeqListTrain
                                if tempListIntraPic(tempIdx) == tempListOrgLib(idxLibPic)
                                    tempPercentInter = 1;
                                else
                                    tempPercentInter = 0;
                                end
                            else
                                tempPercentInter = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(tempIdx, 2);
                            end
                            
                            tempScale = [tempPercentInter tempDiffPic];
                            tempList = [];
                            for idxDqp = 1: numDqp
                                % collect data.
                                tempQPList = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(:, 5);
                                tempUniqueQP = unique(tempQPList);
%                                 if length(tempUniqueQP) < 5
%                                     continue;
%                                 end
                                
                                tempDistL0 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList <= tempUniqueQP(1), 2);
%                                 if length(tempUniqueQP) == 5
%                                     tempDistL1 = tempDistL0;
%                                     tempDistL2 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(2), 2);
%                                     tempDistL3 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(3), 2);
%                                     tempDistL4 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(4), 2);
%                                     tempDistL5 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(5), 2);
%                                 else
%                                     tempDistL1 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(2), 2);
%                                     tempDistL2 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(3), 2);
%                                     tempDistL3 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(4), 2);
%                                     tempDistL4 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(5), 2);
%                                     tempDistL5 = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList == tempUniqueQP(6), 2);
%                                 end

%                                 tempDist(idxDqp, :) = [mean(psnr2mse(tempDistL0)) mean(psnr2mse(tempDistL1)) mean(psnr2mse(tempDistL2)) mean(psnr2mse(tempDistL3)) mean(psnr2mse(tempDistL4)) mean(psnr2mse(tempDistL5))];
                                tempDist(idxDqp) = mean(psnr2mse(tempDistL0));
                                tempList = [tempList bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(:, 2)];
                                tempDistAllLayer(idxDqp) = mean(psnr2mse(bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{idxDqp}{idxLibPic}{idxKeyPic}(tempQPList > tempUniqueQP(1), 2)));
                            end
%                             if length(tempUniqueQP) < 5
%                                 continue;
%                             end
                                
                            eqLinear = @(a,b,x) a.*x+b;
                            tempOpt = fitoptions(eqLinear);
                            tempOpt.StartPoint = [1 1];
                            [tempFitRes, tempGof] = fit(tempDist, tempDistAllLayer, eqLinear, tempOpt);
                            tempScale = [tempScale tempGof.rsquare tempFitRes.a];
%                             for idx = 2: 6
%                                 [tempFitRes, tempGof] = fit(tempDist(:, idx-1), tempDist(:, idx), eqLinear, tempOpt);
%                                 tempScale = [tempScale tempGof.rsquare tempFitRes.a];
%                             end
                            tempScaleAllDQP = [tempScaleAllDQP; tempScale];
                        end
	                end
				end
				tempScaleDAllQP{idxQP} = tempScaleAllDQP;
			end
			scaleD{idxReso} = tempScaleDAllQP;
        end
        
        if flagReadSeqListTrain
            scaleD_Train = scaleD;
        elseif flagReadSeqListTest
            scaleD_Test = scaleD;
        end
	end

	%% test the train model, find the best delta QP and compute the corresponding bd-rate.
	if subFlagTestDQPModelANDComputeBDRate
		subFlagShowResAtAllQPCombination = 1;
		% common models.
		at=0.01;
		listQPOffset = [6 8 8 8];
		% compute rdcost to find the best delta QP.
		% prepare data.
		cellListSeq = cellListSeq_Test;
		numAllSeq = numSeq_Test;
		cellLibInfo = cellLibInfo_Test;
		oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
        bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
        bitPsnrEachLibAllSeq = bitPsnrEachLibAllSeq_Test;
        diffPicAllSeq = diffPicAllSeq_Test;
        cellEncAnalInfo = cellEncAnalInfo_Test;

        % collect data.
        bestDqpAndGainAllLibPicAllQpTrainModel = cell(numReso, 4);
        tempCountAllSeq = 0;
        temp = [];
        tempDqp = [];
        for idxReso = 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
	        % collect data.
	        tempGainAllSeq = zeros(tempNumSeq, 3);
	        tempBestDqpAllSeq = cell(tempNumSeq, 1);
	        tempRdcostAllSeq = cell(tempNumSeq, 1);
            tempVerifyData = cell(numQP, 1);           

	        for idxSeq = 1: tempNumSeq
	        	% prepare data.
	        	tempIP = cellListSeq{idxReso}{idxSeq, 3};
	        	tempSeqLength = cellListSeq{idxReso}{idxSeq, 4};
	        	tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
				tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
	            tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length',tempEstInfoAllLibPic)>0);
				tempEstInfoAllKeyPic = oriContentInfoAllSeq{idxReso}{idxSeq, 1};

				% collect bit psnr.
				tempBitPsnrLibvcAllQP = zeros(numQP, 4);
				tempBitPsnrAnchorAllQP = zeros(numQP, 4);
				tempBestDqpAllQP = zeros(tempNumLibPic, numQP);
	            tempRdcostAllQP = cell(tempNumLibPic, numQP);

				for idxQP = 1: numQP
					% prepare data.
					tempQpKey =listQP(idxQP);
					tempQKey = qp2qstep(tempQpKey);

					% collect data.
					tempSumNumPicLibvc = 0;
					tempSumNumPicAnchor = 0;
                    tempCountKeyPicPerLib = zeros(tempNumLibPic, 1);

					for idxLibPic = 1: tempNumLibPic
						% collect data.
						tempRdcostAllDqp = zeros(numDqp, 1);
						tempSumWDBeta = 0;
                        tempSumPercentInter = [];

						for idxKeyPic = 1: tempNumKeyPic
							% check whether the key picture reference the given library picture.
							idxLibPicForKeyPic = tempListRefLib(idxKeyPic) + 1;
							if idxLibPicForKeyPic ~= idxLibPic
								continue;
							end
							tempKeyIsLib = false;
							% check whether the key picture is the same one as library picture.
							if tempListOrgLib(idxLibPic) == tempListKeyPic(idxKeyPic)
								tempKeyIsLib = true;
							end
							% prepare data.
                            tempCountKeyPicPerLib(idxLibPic) = tempCountKeyPicPerLib(idxLibPic) + 1;
                            tempPercentInterKey = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 2);
                            tempSumPercentInter = [tempSumPercentInter; tempPercentInterKey];

                            if idxKeyPic ~= tempNumKeyPic
                            	tempNumPicPerKey = tempListKeyPic(idxKeyPic+1) - tempListKeyPic(idxKeyPic);
                            else
                            	tempNumPicPerKey = tempSeqLength - tempListKeyPic(idxKeyPic);
                            end
                            
                            if tempKeyIsLib
                                tempBeta = at * tempQKey;
                            else
                                tempBeta = at * tempPercentInterKey * tempQKey;
                            end
                            tempDiffPic = diffPicAllSeq{idxReso}{idxSeq}{idxLibPic}(tempCountKeyPicPerLib(idxLibPic));
                            
                            if tempKeyIsLib
                                %tempW = -0.03*tempDiffPic+0.7267;
                                tempW = -0.0367*tempDiffPic+0.8119;
                            else
                                if tempDiffPic >10
                                    tempW = 0.5;
                                    %tempW = 0;
                                else
                                    %tempW = -0.024*tempDiffPic+1.07;
                                    %tempW = -0.0215*tempDiffPic+0.956;
                                    tempW = 0.95;
                                    %tempW = 0;
                                end
                            end
                            %tempWD = tempW * (tempNumPicPerKey-1) + 1;
                            tempWD = tempW * tempNumPicPerKey;
                            tempSumWDBeta = tempSumWDBeta + tempWD * tempBeta;
                            tempVerifyData{idxQP} = [tempVerifyData{idxQP}; tempPercentInterKey tempDiffPic tempBeta tempW];
                        end
                        if mean(tempSumPercentInter) >0.7
                            tempOffset = -3;
                        else
                            tempOffset = 0;
                        end

                        
                        tempDQP = -3 * log2(tempSumWDBeta) + listQPOffset(idxQP) + tempOffset;
						
                        tempDQP = max(min(round(tempDQP), 0), -11);

                        tempBestDqpIdx = find(listDqp == tempDQP);
	                    tempBestDqpAllQP(idxLibPic, idxQP) = listDqp(tempBestDqpIdx);

						tempBitPsnrLibvc = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{tempBestDqpIdx}{idxLibPic, 1};
						tempNumPicLibvc = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{tempBestDqpIdx}{idxLibPic, 2};
						tempBitPsnrLibvcAllQP(idxQP, :) = tempBitPsnrLibvcAllQP(idxQP, :) + tempBitPsnrLibvc * tempNumPicLibvc;
						tempSumNumPicLibvc = tempSumNumPicLibvc + tempNumPicLibvc;

						tempBitPsnrAnchor = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 1}{idxQP}{idxLibPic, 1};
						tempNumPicAnchor = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 1}{idxQP}{idxLibPic, 2};
						tempBitPsnrAnchorAllQP(idxQP, :) = tempBitPsnrAnchorAllQP(idxQP, :) + tempBitPsnrAnchor * tempNumPicAnchor;
						tempSumNumPicAnchor = tempSumNumPicAnchor + tempNumPicAnchor;
	                    
	                    % collect data.
	                    tempRdcostAllQP{idxLibPic, idxQP} = tempRdcostAllDqp;
	                end
	                %tempCount = tempCount + 1;

					% compute average bitrate and psnr.
					tempBitPsnrLibvcAllQP(idxQP, :) = tempBitPsnrLibvcAllQP(idxQP, :) / tempSumNumPicLibvc;
					tempBitPsnrAnchorAllQP(idxQP, :) = tempBitPsnrAnchorAllQP(idxQP, :) / tempSumNumPicAnchor;
					if tempSumNumPicLibvc ~= tempSumNumPicAnchor
						error('error in computing bit psnr of anchor and libvc');
	                end
	            end
	            tempCountAllSeq = tempCountAllSeq + tempNumLibPic;

				% compute bdrate.
				tempGainAllSeq(idxSeq, :) = bdRateComparation( tempBitPsnrAnchorAllQP, tempBitPsnrLibvcAllQP );
				tempBestDqpAllSeq(idxSeq) = {tempBestDqpAllQP};
                tempDqp = [tempDqp; tempBestDqpAllQP];
	            
	            % collect data.
	            tempRdcostAllSeq{idxSeq} = tempRdcostAllQP;
            end
            temp = [temp; tempGainAllSeq(:, 1)];
	        % collect data.
	        bestDqpAndGainAllLibPicAllQpTrainModel(idxReso, :) = [{tempBestDqpAllSeq} {tempGainAllSeq} {tempRdcostAllSeq} {tempVerifyData}];
        end
        temp

        if subFlagShowResAtAllQPCombination
        	bestDqpAndGainAllLibPicAllQp_BDrate = bestDqpAndGainAllLibPicAllQp_BDrate_Test;
	    	for idxReso = 1: numReso
	    		for idxSeq = 1: length(bestDqpAndGainAllLibPicAllQp_BDrate{idxReso, 3})
	    			for idxLibPic = 1: length(bestDqpAndGainAllLibPicAllQp_BDrate{idxReso, 3}{idxSeq})
	    				for idxCase = 1: length(bestDqpAndGainAllLibPicAllQp_BDrate{idxReso, 3}{idxSeq}{idxLibPic, 2})
	    					plot(bestDqpAndGainAllLibPicAllQp_BDrate{idxReso, 3}{idxSeq}{idxLibPic, 2}{idxCase});
	    				end
	    			end
	    		end
	    	end
	    end
	end

	%% test the train model, find the best delta QP and compute the corresponding bd-rate.
	if subFlagTestTrainModelANDComputeBDRate
		% common models.
		eqDLib = @(a,e,q) a.*e.^0.5.*q;
		eqRLib = @(u,f,b,e,q) u.*exp(-f.*q./(e.^0.5))+b;
		eqDKey = @(a1,a2,e1,e2,e3,s,q1,q2) (1-s).*eqDLib(a2,e1,q1)+s.*eqDLib(a1,e2+eqDLib(a2,e3,q2),q1);
		eqRKey = @(u1,f1,b1,u2,f2,b2,a2,e1,e2,e3,s,q1,q2) (1-s).*eqRLib(u2,f2,b2,e1,q1)+s.*eqRLib(u1,f1,b1,e2+eqDLib(a2,e3,q2),q1);
		% constant coefficients.
		a1=0.1216;a2=0.08302;
		u1=3.578;f1=0.947;b1=0.8733;
		u2=4.645;f2=1.906;b2=0.2732;
		wD=1.647;wR=0.155;

		aDep=[0.6808;0;0.9785;0.7921];
		bDep=[0;0;-0.1386;-0.1401];
        dDep=[-0.15;-0.06;0;0];

		% compute rdcost to find the best delta QP.
		% prepare data.
		cellListSeq = cellListSeq_Test;
		numAllSeq = numSeq_Test;
		cellLibInfo = cellLibInfo_Test;
		oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
        bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
        bitPsnrEachLibAllSeq = bitPsnrEachLibAllSeq_Test;

        % collect data.
        bestDqpAndGainAllLibPicAllQpTrainModel = cell(numReso, 3);
        tempCountAllSeq = 0;
        temp = [];
        for idxReso = 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
	        % collect data.
	        tempGainAllSeq = zeros(tempNumSeq, 3);
	        tempBestDqpAllSeq = cell(tempNumSeq, 1);
	        tempRdcostAllSeq = cell(tempNumSeq, 1);

	        for idxSeq = 1: tempNumSeq
	        	% prepare data.
	        	tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
				tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
	            tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length',tempEstInfoAllLibPic)>0);
				tempEstInfoAllKeyPic = oriContentInfoAllSeq{idxReso}{idxSeq, 1};

				% collect bit psnr.
				tempBitPsnrLibvcAllQP = zeros(numQP, 4);
				tempBitPsnrAnchorAllQP = zeros(numQP, 4);
				tempBestDqpAllQP = zeros(tempNumLibPic, numQP);
	            tempRdcostAllQP = cell(tempNumLibPic, numQP);

				for idxQP = 1: numQP
					% prepare data.
					tempQpKey =listQP(idxQP);
					tempQKey = qp2qstep(tempQpKey);
					tempLambdaPerQP = 0.57 * 2 .^ ((tempQpKey - 12) / 3);

					% collect data.
					tempSumNumPicLibvc = 0;
					tempSumNumPicAnchor = 0;

					for idxLibPic = 1: tempNumLibPic
						% collect data.
						tempRdcostAllDqp = zeros(numDqp, 1);

						for idxKeyPic = 1: tempNumKeyPic
							% check whether the key picture reference the given library picture.
							idxLibPicForKeyPic = tempListRefLib(idxKeyPic) + 1;
							if idxLibPicForKeyPic ~= idxLibPic
								continue;
							end
							tempKeyIsLib = false;
							% check whether the key picture is the same one as library picture.
							if tempListOrgLib(idxLibPic) == tempListKeyPic(idxKeyPic)
								tempKeyIsLib = true;
							end
							% prepare estimated content information.
					        tempEstResiEnergyLib = tempEstInfoAllLibPic{idxLibPic, 1}{1};
	                        if ~tempKeyIsLib
	                            tempEstResiEnergyKeyInter = tempEstInfoAllKeyPic{idxKeyPic, 1}{1};
	                            tempEstResiEnergyKeyIntra = tempEstInfoAllKeyPic{idxKeyPic, 2}{1};
	                            tempPercentInterKey = tempEstInfoAllKeyPic{idxKeyPic, 3};
	                        else
	                            tempEstResiEnergyKeyInter = 0;
	                            tempEstResiEnergyKeyIntra = 0;
	                            tempPercentInterKey = 0;
	                        end

						    % compute rdcost for each delta QP.
						    tempRateCost = 0;
						    tempDistCost = 0;
							for idxDqp = 1: numDqp
								% compute reference QP.
								tempQpLib = tempQpKey + listDqp(idxDqp);
								tempQLib = qp2qstep(tempQpLib);

								% compute distortion.
								tempDistLib = eqDLib(a2,tempEstResiEnergyLib,tempQLib);
								%tempDistKey = eqDKey(a1,a2,tempEstResiEnergyKeyIntra,tempEstResiEnergyKeyInter,tempEstResiEnergyLib,tempPercentInterKey,tempQKey,tempQLib);
								tempS = max(0,aDep(idxReso)*tempPercentInterKey + bDep(idxReso)+dDep(idxQP));
                                tempDistKey = tempDistLib * tempS;
								% compute rate.
								tempRateLib = eqRLib(u2,f2,b2,tempEstResiEnergyLib,tempQLib);
								%tempRateKey = eqRKey(u1,f1,b1,u2,f2,b2,a2,tempEstResiEnergyKeyIntra,tempEstResiEnergyKeyInter,tempEstResiEnergyLib,tempPercentInterKey,tempQKey,tempQLib);
								tempRateKey = 0;
								% compute rdcost.
								if tempKeyIsLib
									%tempRateCost = tempRateLib / wR;
                                    tempRateCost = tempRateLib;
									tempDistCost = tempDistLib;
	                                tempLambda = tempLambdaPerQP;
								else
									tempRateCost = tempRateKey;
									tempDistCost = tempDistKey;
	                                tempLambda = tempLambdaPerQP;
								end
								%tempRdcostAllDqp(idxDqp) = tempRdcostAllDqp(idxDqp) + wD * tempDistCost + tempLambda * wR * tempRateCost;
                                tempRdcostAllDqp(idxDqp) = tempRdcostAllDqp(idxDqp) + tempDistCost + tempLambda * tempRateCost;
	                        end
						end

						% compare rdcost and find the best delta QP.
						[~, tempBestDqpIdx] = min(tempRdcostAllDqp);

						% compute bit psnr with given best delta QP.
	                    %tempBestDqpIdx = find(listDqp == bestDqpAllSeqAllLibPicAllQp_BDrate_Test(idxQP+(tempCountAllSeq+idxLibPic-1)*numQP));
	                    %tempBestDqpAllQP(idxLibPic, idxQP) = listDqp(tempBestDqpIdx);
	                    %tempBestDqpGiven = [-4,-3,-3,-3,-4,-9,-7,-5,-7,-4,-11,-4,-3,-3,-3,-4,-3,-4,-2,-3,-7,-4,-4,-3,-3,-5];
	                    %tempBestDqpGiven = [-4,-3,-3,-3,-4,-9,-7,-5,-7,-6,-11,-6,-4,-4,-4,-6,-6,-6,-2,-6,-10,-6,-7,-5,-5,-8];
	                    tempBestDqpGiven = [-3,0,0,-7,-10,-9,0,-6,-7,-7,-10,-8,-6,-6,-7,-8,-9,0,0,-5,-6,-8,-9,-9,-9,-1;-2,0,0,-6,-9,-9,0,-5,-6,-6,-9,-6,-5,-5,-5,-7,-8,0,0,-3,-5,-6,-7,-8,-8,0;-4,0,-1,-7,-10,-11,0,-6,-7,-7,-11,-8,-5,-6,-6,-8,-9,0,0,-4,-6,-8,-9,-9,-9,-1;-5,-1,-3,-8,-11,-11,-1,-8,-9,-8,-11,-9,-6,-7,-7,-10,-11,-2,0,-5,-8,-9,-11,-11,-11,-3];
                        tempBestDqpIdx = find(listDqp == tempBestDqpGiven(idxQP, tempCountAllSeq+idxLibPic));                   
	                    tempBestDqpAllQP(idxLibPic, idxQP) = listDqp(tempBestDqpIdx);
	                        
	                    
						tempBitPsnrLibvc = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{tempBestDqpIdx}{idxLibPic, 1};
						tempNumPicLibvc = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{tempBestDqpIdx}{idxLibPic, 2};
						tempBitPsnrLibvcAllQP(idxQP, :) = tempBitPsnrLibvcAllQP(idxQP, :) + tempBitPsnrLibvc * tempNumPicLibvc;
						tempSumNumPicLibvc = tempSumNumPicLibvc + tempNumPicLibvc;

						tempBitPsnrAnchor = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 1}{idxQP}{idxLibPic, 1};
						tempNumPicAnchor = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 1}{idxQP}{idxLibPic, 2};
						tempBitPsnrAnchorAllQP(idxQP, :) = tempBitPsnrAnchorAllQP(idxQP, :) + tempBitPsnrAnchor * tempNumPicAnchor;
						tempSumNumPicAnchor = tempSumNumPicAnchor + tempNumPicAnchor;
	                    
	                    % collect data.
	                    tempRdcostAllQP{idxLibPic, idxQP} = tempRdcostAllDqp;
	                end
	                %tempCount = tempCount + 1;

					% compute average bitrate and psnr.
					tempBitPsnrLibvcAllQP(idxQP, :) = tempBitPsnrLibvcAllQP(idxQP, :) / tempSumNumPicLibvc;
					tempBitPsnrAnchorAllQP(idxQP, :) = tempBitPsnrAnchorAllQP(idxQP, :) / tempSumNumPicAnchor;
					if tempSumNumPicLibvc ~= tempSumNumPicAnchor
						error('error in computing bit psnr of anchor and libvc');
	                end
	            end
	            tempCountAllSeq = tempCountAllSeq + tempNumLibPic;

				% compute bdrate.
				tempGainAllSeq(idxSeq, :) = bdRateComparation( tempBitPsnrAnchorAllQP, tempBitPsnrLibvcAllQP );
				tempBestDqpAllSeq(idxSeq) = {tempBestDqpAllQP};
	            
	            % collect data.
	            tempRdcostAllSeq{idxSeq} = tempRdcostAllQP;
            end
            temp = [temp; tempGainAllSeq(:, 1)];
	        % collect data.
	        bestDqpAndGainAllLibPicAllQpTrainModel(idxReso, :) = [{tempBestDqpAllSeq} {tempGainAllSeq} {tempRdcostAllSeq}];
        end
        temp
	end

	%% test the train model, find the best delta QP and compute the corresponding bd-rate.
	if subFlagRdcostOnActualBitDistANDComputeBDRate
		% compute rdcost to find the best delta QP.
		% prepare data.
		cellListSeq = cellListSeq_Test;
		numAllSeq = numSeq_Test;
		cellLibInfo = cellLibInfo_Test;
		oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
        bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
        bitPsnrEachLibAllSeq = bitPsnrEachLibAllSeq_Test;
		bitPsnrAllPicAllLibAllSeq = bitPsnrAllPicAllLibAllSeq_Test;
    
        wD=1.647;wR=0.155;
        % collect data.
        bestDqpAndGainAllLibPicAllQpTrainModel = cell(numReso, 3);
        tempCountAllSeq = 0;
        temp = [];
        for idxReso = 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
	        % collect data.
	        tempGainAllSeq = zeros(tempNumSeq, 3);
	        tempBestDqpAllSeq = cell(tempNumSeq, 1);
	        tempRdcostAllSeq = cell(tempNumSeq, 1);

	        for idxSeq = 1: tempNumSeq
	        	% prepare data.
	        	tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
                tempGOP = cellListSeq{idxReso}{idxSeq, 3};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);

				% collect bit psnr.
				tempBitPsnrLibvcAllQP = zeros(numQP, 4);
				tempBitPsnrAnchorAllQP = zeros(numQP, 4);
				tempBestDqpAllQP = zeros(tempNumLibPic, numQP);
	            tempRdcostAllQP = cell(tempNumLibPic, numQP);

				for idxQP = 1: numQP
					% prepare data.
					tempQpKey =listQP(idxQP);
					tempQKey = qp2qstep(tempQpKey);
					tempLambdaPerQP = 0.57 * 2 .^ ((tempQpKey - 12) / 3);
                    tempLambdaPerQP = 0.4624/16*tempQKey^2;

					% collect data.
					tempSumNumPicLibvc = 0;
					tempSumNumPicAnchor = 0;

					for idxLibPic = 1: tempNumLibPic
						% collect data.
						tempRdcostAllDqp = zeros(numDqp, 1);
						idxKeyPicPerLib = 1;

						for idxKeyPic = 1: tempNumKeyPic
							% check whether the key picture reference the given library picture.
							idxLibPicForKeyPic = tempListRefLib(idxKeyPic) + 1;
							if idxLibPicForKeyPic ~= idxLibPic
								continue;
							end
							tempKeyIsLib = false;
							% check whether the key picture is the same one as library picture.
							if tempListOrgLib(idxLibPic) == tempListKeyPic(idxKeyPic)
								tempKeyIsLib = true;
							end
						    % compute rdcost for each delta QP.
						    tempRateCost = 0;
						    tempDistCost = 0;
							for idxDqp = 1: numDqp
								tempDistKey = psnr2mse(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(idxKeyPicPerLib, 5));
								tempRateKey = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(idxKeyPicPerLib, 4) / tempWidth / tempHeight;
								tempRateLib = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(4) / tempWidth / tempHeight;

								% compute rdcost.
								if tempKeyIsLib
									%tempRateCost = tempRateLib / wR;
                                    tempRateCost = tempRateLib / wR / tempGOP;
									tempDistCost = 0;
	                                tempLambda = tempLambdaPerQP;
								else
									tempRateCost = tempRateKey;
									tempDistCost = tempDistKey;
	                                tempLambda = tempLambdaPerQP;
                                end
                                tempRdcostAllDqp(idxDqp) = tempRdcostAllDqp(idxDqp) + wD * tempDistCost + tempLambda * wR * tempRateCost;
                                %tempRdcostAllDqp(idxDqp) = tempRdcostAllDqp(idxDqp) + tempDistCost + tempLambda * tempRateCost;
	                        end
	                        idxKeyPicPerLib = idxKeyPicPerLib + 1;
						end
						% compute rdcost using actual bit dist of all pictures.
						for idxDqp = 1: numDqp
                            tempRateLib = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(4) / tempWidth / tempHeight;
							tempRateCost = bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic}(:, 1) / tempWidth / tempHeight;
							tempDistCost = psnr2mse(bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic}(:, 2));
                            tempLambdaList = 0.4624/16*bitPsnrAllPicAllLibAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic}(:, 5).^2;
                            tempRdcostAllDqp(idxDqp) = sum(tempDistCost) + sum(tempRateCost .* tempLambdaList) + tempLambdaPerQP * tempRateLib;
						end

						% compare rdcost and find the best delta QP.
						[~, tempBestDqpIdx] = min(tempRdcostAllDqp);

	                    tempBestDqpAllQP(idxLibPic, idxQP) = listDqp(tempBestDqpIdx);
	                        
	                    
						tempBitPsnrLibvc = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{tempBestDqpIdx}{idxLibPic, 1};
						tempNumPicLibvc = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{tempBestDqpIdx}{idxLibPic, 2};
						tempBitPsnrLibvcAllQP(idxQP, :) = tempBitPsnrLibvcAllQP(idxQP, :) + tempBitPsnrLibvc * tempNumPicLibvc;
						tempSumNumPicLibvc = tempSumNumPicLibvc + tempNumPicLibvc;

						tempBitPsnrAnchor = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 1}{idxQP}{idxLibPic, 1};
						tempNumPicAnchor = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 1}{idxQP}{idxLibPic, 2};
						tempBitPsnrAnchorAllQP(idxQP, :) = tempBitPsnrAnchorAllQP(idxQP, :) + tempBitPsnrAnchor * tempNumPicAnchor;
						tempSumNumPicAnchor = tempSumNumPicAnchor + tempNumPicAnchor;
	                    
	                    % collect data.
	                    tempRdcostAllQP{idxLibPic, idxQP} = tempRdcostAllDqp;
	                end
	                %tempCount = tempCount + 1;

					% compute average bitrate and psnr.
					tempBitPsnrLibvcAllQP(idxQP, :) = tempBitPsnrLibvcAllQP(idxQP, :) / tempSumNumPicLibvc;
					tempBitPsnrAnchorAllQP(idxQP, :) = tempBitPsnrAnchorAllQP(idxQP, :) / tempSumNumPicAnchor;
					if tempSumNumPicLibvc ~= tempSumNumPicAnchor
						error('error in computing bit psnr of anchor and libvc');
	                end
	            end
	            tempCountAllSeq = tempCountAllSeq + tempNumLibPic;

				% compute bdrate.
				tempGainAllSeq(idxSeq, :) = bdRateComparation( tempBitPsnrAnchorAllQP, tempBitPsnrLibvcAllQP );
				tempBestDqpAllSeq(idxSeq) = {tempBestDqpAllQP};
	            
	            % collect data.
	            tempRdcostAllSeq{idxSeq} = tempRdcostAllQP;
	        end
	        % collect data.
            temp = [temp; tempGainAllSeq(:, 1)];
	        bestDqpAndGainAllLibPicAllQpTrainModel(idxReso, :) = [{tempBestDqpAllSeq} {tempGainAllSeq} {tempRdcostAllSeq}];
        end
        temp
	end

	%% test the train model, find the best delta QP and compute the corresponding bd-rate.
	if subFlagTestTrainModelOnActualEncInfoANDComputeBDRate
		% common models.
		eqDLib = @(a,e,q) a.*e.^0.5.*q;
		eqRLib = @(u,f,b,e,q) u.*exp(-f.*q./(e.^0.5))+b;
		eqDKey = @(a1,a2,e1,e2,e3,s,q1,q2) (1-s).*eqDLib(a2,e1,q1)+s.*eqDLib(a1,e2+eqDLib(a2,e3,q2),q1);
		eqRKey = @(u1,f1,b1,u2,f2,b2,a2,e1,e2,e3,s,q1,q2) (1-s).*eqRLib(u2,f2,b2,e1,q1)+s.*eqRLib(u1,f1,b1,e2+eqDLib(a2,e3,q2),q1);
		% constant coefficients.
		a1=0.1216;a2=0.08302;
		u1=3.578;f1=0.947;b1=0.8733;
		u2=4.645;f2=1.906;b2=0.2732;
		wD=1.647;wR=0.155;

		aDep=[0.6808;0;0.9785;0.7921];
		bDep=[0;0;-0.1386;-0.1401];
        dDep=[-0.15;-0.06;0;0];
        encResiDep = [0.002 0 0.0025 0.0032]*3;

		% compute rdcost to find the best delta QP.
		% prepare data.
		cellListSeq = cellListSeq_Test;
		numAllSeq = numSeq_Test;
		cellLibInfo = cellLibInfo_Test;
		oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
        bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
        bitPsnrEachLibAllSeq = bitPsnrEachLibAllSeq_Test;

        % collect data.
        bestDqpAndGainAllLibPicAllQpTrainModel = cell(numReso, 3);
        tempCountAllSeq = 0;
        temp = [];
        for idxReso = 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
	        % collect data.
	        tempGainAllSeq = zeros(tempNumSeq, 3);
	        tempBestDqpAllSeq = cell(tempNumSeq, 1);
	        tempRdcostAllSeq = cell(tempNumSeq, 1);

	        for idxSeq = 1: tempNumSeq
	        	% prepare data.
	        	tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
				tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
	            tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length',tempEstInfoAllLibPic)>0);
				tempEstInfoAllKeyPic = oriContentInfoAllSeq{idxReso}{idxSeq, 1};

				% collect bit psnr.
				tempBitPsnrLibvcAllQP = zeros(numQP, 4);
				tempBitPsnrAnchorAllQP = zeros(numQP, 4);
				tempBestDqpAllQP = zeros(tempNumLibPic, numQP);
	            tempRdcostAllQP = cell(tempNumLibPic, numQP);

				for idxQP = 1: numQP
					% prepare data.
					tempQpKey =listQP(idxQP);
					tempQKey = qp2qstep(tempQpKey);
					tempLambdaPerQP = 0.57 * 2 .^ ((tempQpKey - 12) / 3);

					% collect data.
					tempSumNumPicLibvc = 0;
					tempSumNumPicAnchor = 0;

					for idxLibPic = 1: tempNumLibPic
						% collect data.
						tempRdcostAllDqp = zeros(numDqp, 1);

						for idxKeyPic = 1: tempNumKeyPic
							% check whether the key picture reference the given library picture.
							idxLibPicForKeyPic = tempListRefLib(idxKeyPic) + 1;
							if idxLibPicForKeyPic ~= idxLibPic
								continue;
							end
							tempKeyIsLib = false;
							% check whether the key picture is the same one as library picture.
							if tempListOrgLib(idxLibPic) == tempListKeyPic(idxKeyPic)
								tempKeyIsLib = true;
							end
							% prepare estimated content information.
					        tempEstResiEnergyLib = tempEstInfoAllLibPic{idxLibPic, 1}{1};
	                        if ~tempKeyIsLib
	                            tempEstResiEnergyKeyInter = tempEstInfoAllKeyPic{idxKeyPic, 1}{1};
	                            tempEstResiEnergyKeyIntra = tempEstInfoAllKeyPic{idxKeyPic, 2}{1};
	                            %tempPercentInterKey = tempEstInfoAllKeyPic{idxKeyPic, 3};
	                        else
	                            tempEstResiEnergyKeyInter = 0;
	                            tempEstResiEnergyKeyIntra = 0;
	                            tempPercentInterKey = 0;
	                        end

						    % compute rdcost for each delta QP.
						    tempRateCost = 0;
						    tempDistCost = 0;
							for idxDqp = 1: numDqp
								% compute reference QP.
								tempQpLib = tempQpKey + listDqp(idxDqp);
								tempQLib = qp2qstep(tempQpLib);

								% compute distortion.
								tempDistLib = eqDLib(a2,tempEstResiEnergyLib,tempQLib);
								tempRateLib = eqRLib(u2,f2,b2,tempEstResiEnergyLib,tempQLib);
								tempRateKey = 0;
	                            tempPercentInterKey = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 2);
	                            tempDistKey = tempPercentInterKey * encResiDep(idxReso) * tempQKey * 0.9 * tempDistLib;
								% compute rdcost.
								if tempKeyIsLib
									%tempRateCost = tempRateLib / wR;
                                    tempRateCost = tempRateLib;
									tempDistCost = 0;
	                                tempLambda = tempLambdaPerQP;
								else
									tempRateCost = tempRateKey;
									tempDistCost = tempDistKey;
	                                tempLambda = tempLambdaPerQP;
								end
								tempRdcostAllDqp(idxDqp) = tempRdcostAllDqp(idxDqp) + wD * tempDistCost + tempLambda * wR * tempRateCost;
                                %tempRdcostAllDqp(idxDqp) = tempRdcostAllDqp(idxDqp) + tempDistCost + tempLambda * tempRateCost;
	                        end
						end

						% compare rdcost and find the best delta QP.
						[~, tempBestDqpIdx] = min(tempRdcostAllDqp);

						% compute bit psnr with given best delta QP.
	                    %tempBestDqpIdx = find(listDqp == bestDqpAllSeqAllLibPicAllQp_BDrate_Test(idxQP+(tempCountAllSeq+idxLibPic-1)*numQP));
	                    %tempBestDqpAllQP(idxLibPic, idxQP) = listDqp(tempBestDqpIdx);
	                    %tempBestDqpGiven = [-4,-3,-3,-3,-4,-9,-7,-5,-7,-4,-11,-4,-3,-3,-3,-4,-3,-4,-2,-3,-7,-4,-4,-3,-3,-5];
	                    %tempBestDqpGiven = [-4,-3,-3,-3,-4,-9,-7,-5,-7,-6,-11,-6,-4,-4,-4,-6,-6,-6,-2,-6,-10,-6,-7,-5,-5,-8];
	                    %tempBestDqpGiven = [-3,-3,-3,-3,-5,-10,-5,-4,-5,-4,-11,-5,-3,-4,-4,-6,-6,-5,-2,-3,-6,-5,-6,-5,-4,-5];
                        %tempBestDqpIdx = find(listDqp == tempBestDqpGiven(tempCountAllSeq+idxLibPic));                      
	                    tempBestDqpAllQP(idxLibPic, idxQP) = listDqp(tempBestDqpIdx);
	                        
	                    
						tempBitPsnrLibvc = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{tempBestDqpIdx}{idxLibPic, 1};
						tempNumPicLibvc = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 2}{idxQP}{tempBestDqpIdx}{idxLibPic, 2};
						tempBitPsnrLibvcAllQP(idxQP, :) = tempBitPsnrLibvcAllQP(idxQP, :) + tempBitPsnrLibvc * tempNumPicLibvc;
						tempSumNumPicLibvc = tempSumNumPicLibvc + tempNumPicLibvc;

						tempBitPsnrAnchor = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 1}{idxQP}{idxLibPic, 1};
						tempNumPicAnchor = bitPsnrEachLibAllSeq{idxReso}{idxSeq, 1}{idxQP}{idxLibPic, 2};
						tempBitPsnrAnchorAllQP(idxQP, :) = tempBitPsnrAnchorAllQP(idxQP, :) + tempBitPsnrAnchor * tempNumPicAnchor;
						tempSumNumPicAnchor = tempSumNumPicAnchor + tempNumPicAnchor;
	                    
	                    % collect data.
	                    tempRdcostAllQP{idxLibPic, idxQP} = tempRdcostAllDqp;
	                end
	                %tempCount = tempCount + 1;

					% compute average bitrate and psnr.
					tempBitPsnrLibvcAllQP(idxQP, :) = tempBitPsnrLibvcAllQP(idxQP, :) / tempSumNumPicLibvc;
					tempBitPsnrAnchorAllQP(idxQP, :) = tempBitPsnrAnchorAllQP(idxQP, :) / tempSumNumPicAnchor;
					if tempSumNumPicLibvc ~= tempSumNumPicAnchor
						error('error in computing bit psnr of anchor and libvc');
	                end
	            end
	            tempCountAllSeq = tempCountAllSeq + tempNumLibPic;

				% compute bdrate.
				tempGainAllSeq(idxSeq, :) = bdRateComparation( tempBitPsnrAnchorAllQP, tempBitPsnrLibvcAllQP );
				tempBestDqpAllSeq(idxSeq) = {tempBestDqpAllQP};
	            
	            % collect data.
	            tempRdcostAllSeq{idxSeq} = tempRdcostAllQP;
	        end
	        % collect data.
            temp = [temp; tempGainAllSeq(:, 1)];
	        bestDqpAndGainAllLibPicAllQpTrainModel(idxReso, :) = [{tempBestDqpAllSeq} {tempGainAllSeq} {tempRdcostAllSeq}];
        end
        temp
	end

	if subFlagCheckEncInfoVSEstiInfo
		% compute rdcost to find the best delta QP.
		% prepare data.
		cellListSeq = cellListSeq_Test;
		numAllSeq = numSeq_Test;
		cellLibInfo = cellLibInfo_Test;
		cellEncAnalInfo = cellEncAnalInfo_Test;
		oriContentInfoAllSeq = oriContentInfoAllSeq_Test;
        bitPsnrLibKeyAllSeq = bitPsnrLibKeyAllSeq_Test;
        bitPsnrEachLibAllSeq = bitPsnrEachLibAllSeq_Test;

        % collect data.
        encAnalInfoVSEstiInfo = cell(numReso, 1);
        encPercent = [];
        encPercentSum = [];
        encResiVSEstiResi = [];
        fitResResiDep = cell(numReso, 1);
        for idxReso = 1: numReso
        	if numAllSeq(idxReso) < 1
				continue;
			end
			% prepare data.
			tempNumSeq = numAllSeq(idxReso);
	        % collect data.
	        tempEncAnalInfoVSEstiInfo = cell(tempNumSeq, 1);
	        tempFitResResiDepAllSeq = cell(tempNumSeq, 1);

	        for idxSeq = 1: tempNumSeq
	        	% prepare data.
	        	tempWidth = cellListSeq{idxReso}{idxSeq, 5};
	            tempHeight = cellListSeq{idxReso}{idxSeq, 6};
				tempListKeyPic = cellLibInfo{idxReso}{idxSeq, 4};
				tempListOrgLib = cellLibInfo{idxReso}{idxSeq, 5};
				tempListRefLib = cellLibInfo{idxReso}{idxSeq, 6};
				tempNumLibPic = length(tempListOrgLib);
				tempNumKeyPic = length(tempListKeyPic);
				tempEstInfoAllLibPic = oriContentInfoAllSeq{idxReso}{idxSeq, 2};
	            tempEstInfoAllLibPic = tempEstInfoAllLibPic(cellfun('length',tempEstInfoAllLibPic)>0);
				tempEstInfoAllKeyPic = oriContentInfoAllSeq{idxReso}{idxSeq, 1};

				% collect bit psnr.
				tempEncAnalInfoVSEstiInfoAllLib = cell(tempNumLibPic, 2);
				tempFitResResiDepAllLib = cell(tempNumLibPic, 1);

				for idxLibPic = 1: tempNumLibPic
					% collect library data.
					tempEncAnalInfoVSEstiInfoPerLib = [];
					tempEstLibInfo = tempEstInfoAllLibPic{idxLibPic, 1}{1};
					for idxQP = 1: numQP
						for idxDqp = 1: numDqp
							tempEncAnalInfo = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{1}(idxLibPic, 3);
							tempActualDist = psnr2mse(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(5));
							tempActualRate = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(4) / tempWidth / tempHeight;
							tempQ = qp2qstep(listQP(idxQP) + listDqp(idxDqp));
							tempEstiDist = tempEncAnalInfo*eqT(2/3, tempQ/((tempEncAnalInfo/2)^0.5));
							tempEstiRate = equationEntropy(2/3, tempQ/((tempEncAnalInfo/2)^0.5));
							tempEncAnalInfoVSEstiInfoPerLib = [tempEncAnalInfoVSEstiInfoPerLib; tempEstLibInfo tempEncAnalInfo tempEstiDist tempActualDist tempEstiRate tempActualRate];
						end
					end
					tempEncAnalInfoVSEstiInfoPerLib = unique(tempEncAnalInfoVSEstiInfoPerLib, 'rows');
					tempEncAnalInfoVSEstiInfoAllLib{idxLibPic, 1} = tempEncAnalInfoVSEstiInfoPerLib;
                    
                    encResiVSEstiResi = [encResiVSEstiResi; cellEncAnalInfo{idxReso}{idxSeq}{4}{1}{1}(idxLibPic, 3) tempEstLibInfo];
                    
                    % collect key data.
                    tempNumKeyPicPerLib = sum(tempListRefLib == idxLibPic - 1);
                    tempEncAnalInfoVSEstiInfoPerLibAllKey = cell(tempNumKeyPicPerLib - 1, numQP);
					tempKeyCountPerLib = 1;
                    tempKeyCountPerLib2 = 1;
                    tempEncPercentSum = 0;
                    tempFitResAllKeyPic = [];
			        for idxKeyPic = 1: tempNumKeyPic
			        	if tempListRefLib(idxKeyPic) ~= idxLibPic - 1
			        		continue;
			        	end
			        	if tempListKeyPic(idxKeyPic) == tempListOrgLib(idxLibPic)
                            tempKeyCountPerLib2 = tempKeyCountPerLib2 + 1;
			        		continue;
			        	end
			        	tempEstResiEnergyKeyInter = tempEstInfoAllKeyPic{idxKeyPic, 1}{1};
				        tempEstResiEnergyKeyIntra = tempEstInfoAllKeyPic{idxKeyPic, 2}{1};
				        tempPercentInterKey = tempEstInfoAllKeyPic{idxKeyPic, 3};

				        % collect train data.
				        tempDataAllQP = [];

						for idxQP = 1: numQP
							% collect data.
				        	tempEncAnalInfoVSEstiInfoPerLibPerKey = [];
				        	tempDataAllDqp = [];
							for idxDqp = 1: numDqp
                                tempLibDist = psnr2mse(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 1}(5));
                                tempActualKeyDist = psnr2mse(bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(tempKeyCountPerLib2, 5));
                                tempActualKeyRate = bitPsnrLibKeyAllSeq{idxReso}{idxSeq}{idxQP}{idxDqp}{idxLibPic, 2}(tempKeyCountPerLib2, 4) / tempWidth / tempHeight;
								tempEncAnalInfoPercent = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 2);
                                tempEncInterCount = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 4);
                                tempEncIntraCount = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 7);
                                if tempEncInterCount ~= 0
                                    tempEncAnalInfoInterResiVar = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 5) * tempWidth * tempHeight / tempEncInterCount;
                                else
                                    tempEncAnalInfoInterResiVar = 0;
                                end
                                if tempEncIntraCount ~= 0
                                    tempEncAnalInfoIntraResiVar = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{idxDqp}{2}(idxKeyPic, 8) * tempWidth * tempHeight / tempEncIntraCount;
                                else
                                	tempEncAnalInfoIntraResiVar = 0;
                                end
								tempEncAnalInfoVSEstiInfoPerLibPerKey = [tempEncAnalInfoVSEstiInfoPerLibPerKey; tempPercentInterKey tempEncAnalInfoPercent tempEstResiEnergyKeyInter tempLibDist tempEncAnalInfoInterResiVar tempActualKeyDist tempActualKeyRate tempEstResiEnergyKeyIntra tempEncAnalInfoIntraResiVar];
								if idxQP == numQP
                                    tempDataAllDqp = [tempDataAllDqp; tempLibDist tempActualKeyDist tempEncAnalInfoInterResiVar];
                                end
                            end

                            % collect original residue.
                            tempEncAnalInfoPercentQkQ0 = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 2);
                            tempEncInterCount = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 4);
                            tempEncIntraCount = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 7);
                            if tempEncInterCount ~= 0
                                tempEncAnalInfoInterResiVarQkQ0 = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 5) * tempWidth * tempHeight / tempEncInterCount;
                            else 
                                tempEncAnalInfoInterResiVarQkQ0 = 0;
                            end
                            if tempEncIntraCount ~= 0
                                tempEncAnalInfoIntraResiVarQkQ0 = cellEncAnalInfo{idxReso}{idxSeq}{idxQP}{numDqp+1}{2}(idxKeyPic, 8) * tempWidth * tempHeight / tempEncIntraCount;
                            else
                                tempEncAnalInfoIntraResiVarQkQ0 = 0;
                            end
                            tempEncAnalInfoVSEstiInfoPerLibPerKey = [tempEncAnalInfoVSEstiInfoPerLibPerKey; tempPercentInterKey tempEncAnalInfoPercentQkQ0 tempEstResiEnergyKeyInter 0 tempEncAnalInfoInterResiVarQkQ0 0 0 tempEstResiEnergyKeyIntra tempEncAnalInfoIntraResiVarQkQ0];
							
							% collect train data.
                            tempDataAllDqp = [tempDataAllDqp tempEncAnalInfoInterResiVarQkQ0*ones(size(tempDataAllDqp,1),1) tempEncAnalInfoPercentQkQ0*ones(size(tempDataAllDqp,1),1)];
                            tempDataAllQP = [tempDataAllQP; tempDataAllDqp];
							
							% collect original residue on original pixel.
                            tempEncAnalInfoPercentQ0Q0 = cellEncAnalInfo{idxReso}{idxSeq}{numQP+1}{2}(idxKeyPic, 2);
                            tempEncInterCount = cellEncAnalInfo{idxReso}{idxSeq}{numQP+1}{2}(idxKeyPic, 4);
                            tempEncIntraCount = cellEncAnalInfo{idxReso}{idxSeq}{numQP+1}{2}(idxKeyPic, 7);
                            if tempEncInterCount ~= 0
                                tempEncAnalInfoInterResiVarQ0Q0 = cellEncAnalInfo{idxReso}{idxSeq}{numQP+1}{2}(idxKeyPic, 5) * tempWidth * tempHeight / tempEncInterCount;
                            else
                                tempEncAnalInfoInterResiVarQ0Q0 = 0;
                            end
                            if tempEncIntraCount ~= 0
                                tempEncAnalInfoIntraResiVarQ0Q0 = cellEncAnalInfo{idxReso}{idxSeq}{numQP+1}{2}(idxKeyPic, 8) * tempWidth * tempHeight / tempEncIntraCount;
                            else
                                tempEncAnalInfoIntraResiVarQ0Q0 = 0;
                            end
                            tempEncAnalInfoVSEstiInfoPerLibPerKey = [tempEncAnalInfoVSEstiInfoPerLibPerKey; tempPercentInterKey tempEncAnalInfoPercentQ0Q0 tempEstResiEnergyKeyInter 0 tempEncAnalInfoInterResiVarQ0Q0 0 0 tempEstResiEnergyKeyIntra tempEncAnalInfoIntraResiVarQ0Q0];
                            
                            tempEncAnalInfoVSEstiInfoPerLibAllKey{tempKeyCountPerLib, idxQP} = tempEncAnalInfoVSEstiInfoPerLibPerKey;
                        end
                        encPercent = [encPercent;tempEncAnalInfoPercentQkQ0 tempEncAnalInfoInterResiVarQkQ0];
                        if tempEncAnalInfoInterResiVarQkQ0 ~= 0
                            tempEncPercentSum = tempEncPercentSum + tempEncAnalInfoPercentQkQ0;
                        end
                        
						tempKeyCountPerLib = tempKeyCountPerLib + 1;
                        tempKeyCountPerLib2 = tempKeyCountPerLib2 + 1;

                        % train data.
                        eqResiDep = @(a,b,x) a.*x+b;
                        tempXList = tempDataAllQP(:, 1) .* tempDataAllQP(:, 4) .* tempDataAllQP(:, 5);
                        tempYList = tempDataAllQP(:, 2);
                        [tempXList, tempYList] = prepareCurveData( tempXList, tempYList );
                        if isempty(tempXList)
                            continue;
                        end
                        [tempFitRes, tempGof] = fit(tempXList, tempYList, eqResiDep);
                        %plot(tempFitRes, tempXList, tempYList);
                        tempFitResAllKeyPic = [tempFitResAllKeyPic; tempGof.rsquare tempFitRes.a tempFitRes.b];
                    end
                    encPercentSum = [encPercentSum; tempEncPercentSum];
					tempEncAnalInfoVSEstiInfoAllLib{idxLibPic, 2} = tempEncAnalInfoVSEstiInfoPerLibAllKey;
					tempFitResResiDepAllLib{idxLibPic} = tempFitResAllKeyPic;
				end
	            % collect data.
	            tempEncAnalInfoVSEstiInfo{idxSeq} = tempEncAnalInfoVSEstiInfoAllLib;
	            tempFitResResiDepAllSeq{idxSeq} = tempFitResResiDepAllLib;
	        end
	        % collect data.
	        encAnalInfoVSEstiInfo{idxReso} = tempEncAnalInfoVSEstiInfo;
	        fitResResiDep{idxReso} = tempFitResResiDepAllSeq;
	    end
	end
end

if flagPlotIdealCurve
    close all;
    % show laplacian distribution.
    equationLaplacian = @(c,x) ((2./c).^0.5./2.*exp(-(2./c).^0.5.*abs(x)));
    coeffList = -10:1:10;
    resiEner = [1 5 15 30 50]';
    figure;
    for idx = 1: length(resiEner)
        plot(coeffList, equationLaplacian(resiEner(idx), coeffList));
        hold on;
    end
    title('laplacian distribution');
    xlabel('coeff'); ylabel('probability');
    hold off;

    % show integral.
    equationIntegral = @(a,b,x) (-0.5.*exp(-a.*x).*((a.*x+1).^2+1)+1+1./(1-exp(-x)).*(-0.5*exp(-(1+a).*x)).*(((1-b).*x+1).^2+1-exp(x).*((b.*x-1).^2+1)));
    cqList = [0:0.1:1 2:1:20];
    cqList = [4:1:20];
    aIntra = 2/3; bIntra = 1-aIntra;
    aInter = 5/6; bInter = 1-aInter;
    figure;
    plot(cqList, equationIntegral(aIntra, bIntra, cqList), '-r','LineWidth',2);
    hold on;
    plot(cqList, equationIntegral(aInter, bInter, cqList), '-','LineWidth',2);
    title('T(x)');
    xlabel('Q/((resiEnergy/2)^0.5)'); ylabel('normalized distortion(D/resiEnergy)');
    legend('intra','inter');
    hold off;

    % show entropy.
    equationEntropy = @(a,x) (exp(-a.*x).*(log2((1-exp(-a.*x))./(1-exp(-x)))+1+x.*(a.*(1-exp(-x))+exp(-x))./(1-exp(-x))/log(2))-log2(1-exp(-a.*x)));
    equationEntropySimp = @(d,x) (1+d).*exp(-d.*x);
    figure;
    plot(cqList, equationEntropy(aIntra, cqList), '-r','LineWidth',2);
    hold on;
    plot(cqList, equationEntropy(aInter, cqList), '-', 'LineWidth',2);
    plot(cqList, equationEntropySimp(aIntra, cqList), 'd', 'LineWidth',2);
    plot(cqList, equationEntropySimp(aInter, cqList), 'd', 'LineWidth',2);
    title('H(x)');
    xlabel('Q/((resiEnergy/2)^0.5)'); ylabel('bpp(bits per pixel)');
    legend('intra','inter');
    hold off;
end