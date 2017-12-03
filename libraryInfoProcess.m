%% This script is used to read library information for all sequences with all QPs and delta QPs.
flagReadSeqListTest = 0;
flagReadSeqListTrain = 0;
flagReadLibInfo = 0;
flagReadLogFile = 0;
flagPrepareContentInfo = 0;
flagFitForBestDqpEachLib = 0;
flagCheckResultValidity = 0;
flagModelVerify = 1;

% prepare basic information.
listQP = 22:5:37;
listDqp = 2:-1:-11;
numQP = length(listQP);
numDqp = length(listDqp);

%% read sequences list file of CTC.
if flagReadSeqListTest
	% display information.
	disp('reading sequence list of CTC.\n');

    rcVecSeqsB = readSeqList('.\\SeqList_B.txt');
    rcVecSeqsC = readSeqList('.\\SeqList_C.txt');
    rcVecSeqsD = readSeqList('.\\SeqList_D.txt');
    numVecSeqsB = size(rcVecSeqsB, 1);
    numVecSeqsC = size(rcVecSeqsC, 1);
    numVecSeqsD = size(rcVecSeqsD, 1);

    % collect data.
    cellListSeq_Test = [rcVecSeqsB; rcVecSeqsC ; rcVecSeqsD];
    numAllSeq_Test = numVecSeqsB + numVecSeqsC + numVecSeqsD;

    % set path.
    path = '.\log-ctc\';
    pathYuv = '\\mcl\MCL_Space\TestSeq\';
    cellListSeq = cellListSeq_Test;
    numAllSeq = numAllSeq_Test;
    
    % clear temp data.
    clear rcVecSeqsB rcVecSeqsC rcVecSeqsD numVecSeqsB numVecSeqsC numVecSeqsD;
end

%% read sequences list file of train sequences.
if flagReadSeqListTrain
    % display information.
	disp('reading train sequence list.\n');

    rcVecSeqs = readSeqList('.\\SeqList-train-all.txt');
    numVecSeqs = size(rcVecSeqs, 1);

    % collect data.
    cellListSeq_Train = rcVecSeqs;
    numAllSeq_Train = numVecSeqs;
    
    % set path.
    path = '.\log-train\';
    pathYuv = '\\mcl\MCL_Space\TestSeq\';
    cellListSeq = cellListSeq_Train;
    numAllSeq = numAllSeq_Train;

    % clear temp data.
    clear rcVecSeqs numVecSeqs;
end

%% read library information file.
if flagReadLibInfo
	% display information.
	disp('reading library information of sequences.\n');
	% collect data.
	cellLibInfo = cell(numAllSeq, 7);
	for idxSeq = 1: numAllSeq
		% for one sequence, libInfo is same for all qp and dqp, except the area of inter prediction.
		fileLibInfo = [path 'analysis\libraryInfo_' cellListSeq{idxSeq, 2} '.txt'];
		libInfo = readLibInfo( fileLibInfo );

		listOrgLib = cell2mat(libInfo(4));
		listIntraPic = cell2mat(libInfo(11));
		listRefLib = cell2mat(libInfo(12));
		listLibFreq = cell2mat(libInfo(5));
		listCenterCost = cell2mat(libInfo(7));
		listClusterCost = cell2mat(libInfo(8));

		cellLibInfo(idxSeq, 1) = {listCenterCost};
		cellLibInfo(idxSeq, 2) = {listClusterCost};
		cellLibInfo(idxSeq, 3) = {listLibFreq};
		cellLibInfo(idxSeq, 4) = {listIntraPic};
		cellLibInfo(idxSeq, 5) = {listOrgLib};
		cellLibInfo(idxSeq, 6) = {listRefLib};
	end % seq

	% collect data for train or test.
	if flagReadSeqListTrain
		cellLibInfo_Train = cellLibInfo;
	elseif flagReadSeqListTest
		cellLibInfo_Test = cellLibInfo;
	end

	% clear temp data.
	clear dqp qp strDqp strQP fileLibInfo libInfo listSCD listOrgLib listIntraPic listRefLib listLibFreq listCenterCost listClusterCost;
end

%% read log file.
if flagReadLogFile
	subFlagReadSeqLog = 0;
	subFlagReadKeyLog = 0;
	subFlagReadBitPsnrEachLib = 1;
	subFlagReadRDcostEachLib = 1;

	%% read sequence encoding log file for bitrate and psnr.
	if subFlagReadSeqLog
		gainYUVAllSeq = zeros(numAllSeq, numDqp);
		for idxSeq = 1: numAllSeq
			% read anchor.
			bitPsnrAnchor = zeros(4, 5);
			for idxQP = 1: numQP
				strQP = num2str(listQP(idxQP));
				fileAnchor = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_anchor_gop16_enc.txt'];
				bitPsnrAnchor(idxQP, :) = readSumBitPsnr(fileAnchor);
			end

			% read libvc.
			for idxDqp = 1: numDqp
				if(listDqp(idxDqp) <= 0)
					strDqp = num2str(-listDqp(idxDqp));
				else
					strDqp = ['0' num2str(listDqp(idxDqp))];
				end
				bitPsnrLib = zeros(4, 5);
				for idxQP = 1: numQP
					strQP = num2str(listQP(idxQP));
					fileLib = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_libvc_gop16_fixQP_' strDqp '_enc.txt'];
					bitPsnrLib(idxQP, :) = readSumBitPsnr(fileLib);
				end
				gain = bdRateComparation( bitPsnrAnchor, bitPsnrLib );
				gainYUVAllSeq(idxSeq, idxDqp) = gain(4);
			end
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			gainYUVAllSeq_Train = gainYUVAllSeq;
		elseif flagReadSeqListTest
			gainYUVAllSeq_Test = gainYUVAllSeq;
		end

		% clear temp data.
		clear bitPsnrAnchor strQP fileAnchor strDqp bitPsnrLib fileLib gain;
	end

	%% read sequence encoding log file for BD-rate of key frames.
	if subFlagReadKeyLog
		gainKeyYUVAllSeq = zeros(numAllSeq, numDqp);
		for idxSeq = 1: numAllSeq
			% read anchor.
			bitPsnrAnchor = zeros(4, 5);
			for idxQP = 1: numQP
				strQP = num2str(listQP(idxQP));
				fileAnchor = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_anchor_gop16_enc.txt'];
				bitPsnrAnchor(idxQP, :) = readKeyFrameSumBitPsnr(fileAnchor);
			end

			for idxDqp = 1: numDqp
				if(listDqp(idxDqp) <= 0)
					strDqp = num2str(-listDqp(idxDqp));
				else
					strDqp = ['0' num2str(listDqp(idxDqp))];
				end
				bitPsnrLib = zeros(4, 5);
				for idxQP = 1: numQP
					strQP = num2str(listQP(idxQP));
					fileLib = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_libvc_gop16_fixQP_' strDqp '_enc.txt'];
					bitPsnrLib(idxQP, :) = readKeyFrameSumBitPsnr(fileLib);
				end
				gain = bdRateComparation( bitPsnrAnchor, bitPsnrLib );
				gainKeyYUVAllSeq(idxSeq, idxDqp) = gain(4);
			end
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			gainKeyYUVAllSeq_Train = gainKeyYUVAllSeq;
		elseif flagReadSeqListTest
			gainKeyYUVAllSeq_Test = gainKeyYUVAllSeq;
		end

		% clear temp data.
		clear bitPsnrAnchor strQP fileAnchor strDqp bitPsnrLib fileLib gain;
	end

	%% read average bitrate and psnr of frames for each library picture.
	if subFlagReadBitPsnrEachLib
		% display information.
		disp('reading bitrate and PSNR for each cluster of one library picture.\n');
		% collect data.
		gainYEachLibAllSeq = cell(numAllSeq, numDqp);
		gainYEachLibAllSeqKeyFrames = cell(numAllSeq, numDqp);
	    numFrameEachLibAllSeq = cell(numAllSeq, 1);
		for idxSeq = 1: numAllSeq
	        frameRate = cellListSeq{idxSeq, 3};
	    	% prepare data.
			listIntraPic = cellLibInfo{idxSeq, 4};
			listOrgLib = cellLibInfo{idxSeq, 5};
			listRefLib = cellLibInfo{idxSeq, 6};
	    	numLibPic = length(listOrgLib);
			% read anchor.
			% collect data.
			bitPsnrAnchorAllLib = zeros(4* numLibPic, 5);
			bitPsnrAnchorAllLibKeyFrames = zeros(4* numLibPic, 4);
			for idxQP = 1: numQP
				strQP = num2str(listQP(idxQP));
				fileAnchor = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_anchor_gop16_enc.txt'];
				[bitPsnrAnchor, bitPsnrAnchorKeyFrames, ~, ~, listNumFrames] = readBitPsnrSingleLib(fileAnchor, listOrgLib, listIntraPic, listRefLib, frameRate, 0);
				bitPsnrAnchorAllLib(idxQP: 4: end, :) = bitPsnrAnchor;
				bitPsnrAnchorAllLibKeyFrames(idxQP: 4: end, :) = bitPsnrAnchorKeyFrames;
			end
			numFrameEachLibAllSeq(idxSeq) = {listNumFrames};

			for idxDqp = 1: numDqp
				% prepare path.
				if(listDqp(idxDqp) <= 0)
					strDqp = num2str(-listDqp(idxDqp));
				else
					strDqp = ['0' num2str(listDqp(idxDqp))];
				end
				% collect data.
				bitPsnrLibAllLib = zeros(4* numLibPic, 5);
				bitPsnrLibAllLibKeyFrames = zeros(4* numLibPic, 4);
				for idxQP = 1: numQP
					strQP = num2str(listQP(idxQP));
					fileLib = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_libvc_gop16_fixQP_' strDqp '_enc.txt'];
					[bitPsnrLib, bitPsnrLibKeyFrames, ~, ~, listNumFrames] = readBitPsnrSingleLib(fileLib, listOrgLib, listIntraPic, listRefLib, frameRate, 1);
					bitPsnrLibAllLib(idxQP: 4: end, :) = bitPsnrLib;
					bitPsnrLibAllLibKeyFrames(idxQP: 4: end, :) = bitPsnrLibKeyFrames;
				end

				% compute BD-rate for each group of library picture.
				gain = zeros(numLibPic, 4);
				gainKeyFrames = zeros(numLibPic, 3);
				for idxLibPic = 1: numLibPic
					% BD-rate for all frames.
					bitPsnrAnchor = bitPsnrAnchorAllLib(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
					bitPsnrLib = bitPsnrLibAllLib(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
					gain(idxLibPic, :) = bdRateComparation( bitPsnrAnchor, bitPsnrLib );
					% BD-rate for key frames.
					bitPsnrAnchorKeyFrames = bitPsnrAnchorAllLibKeyFrames(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
					bitPsnrLibKeyFrames = bitPsnrLibAllLibKeyFrames(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
					gainKeyFrames(idxLibPic, :) = bdRateComparation( bitPsnrAnchorKeyFrames, bitPsnrLibKeyFrames );
				end
				gainYEachLibAllSeq(idxSeq, idxDqp) = {gain(:, 1)};
				gainYEachLibAllSeqKeyFrames(idxSeq, idxDqp) = {gainKeyFrames(:, 1)};
			end
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			gainYEachLibAllSeq_Train = gainYEachLibAllSeq;
			gainYEachLibAllSeqKeyFrames_Train = gainYEachLibAllSeqKeyFrames;
			numFrameEachLibAllSeq_Train = numFrameEachLibAllSeq;
		elseif flagReadSeqListTest
			gainYEachLibAllSeq_Test = gainYEachLibAllSeq;
			gainYEachLibAllSeqKeyFrames_Test = gainYEachLibAllSeqKeyFrames;
			numFrameEachLibAllSeq_Test = numFrameEachLibAllSeq;
		end

		% clear temp data.
		clear frameRate dqp qp strDqp strQP fileLibInfo numLibPic strDqp fileLib;
		clear strQP fileAnchor listLibPoc listNumFrames;
		clear gain bitPsnrAnchor bitPsnrLib bitPsnrAnchorKeyFrames bitPsnrLibKeyFrames;
		clear bitPsnrAnchorAllLib bitPsnrLibAllLib bitPsnrAnchorAllLibKeyFrames bitPsnrLibAllLibKeyFrames;
	end

	%% read RDcost for each cluster containning single library picture.
	if subFlagReadRDcostEachLib
		% display information.
		disp('reading RDcost for each cluster of one library picture.\n');

		% prepare basic information.
		lambdaFactor = 0.57;

		% collect data for each sequence, each library picture, each base qp, each dqp.
		cellRDcostAllSeq = cell(numAllSeq, 3); % 1 for anchor , 2 for libvc, 3 for numframes of each library picture.
		cellRDcostAllSeqKeyFrames = cell(numAllSeq, 3);

		% compute RDcost.
		for idxSeq = 1: numAllSeq
			% prepare data.
			listLibFreq = cellLibInfo{idxSeq, 3};
			listIntraPic = cellLibInfo{idxSeq, 4};
			listOrgLib = cellLibInfo{idxSeq, 5};
			listRefLib = cellLibInfo{idxSeq, 6};

			numLibPic = length(listOrgLib);
			widthSeq = cellListSeq{idxSeq, 5};
			heightSeq = cellListSeq{idxSeq, 6};

			frameRate = cellListSeq{idxSeq, 3};

			% read anchor.
			% collect data.
			cellRDcostAnchorAllLib = cell(numLibPic, numQP);
			cellRDcostAnchorAllLibKeyFrames = cell(numLibPic, numQP);

			for idxQP = 1: numQP
				% collect data.
				strQP = num2str(listQP(idxQP));
				fileAnchor = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_anchor_gop16_enc.txt'];
				[bitPsnrAnchor, bitPsnrAnchorKeyFrames, sseAnchor, sseAnchorKeyFrames, listNumFrames] = readBitPsnrSingleLib(fileAnchor, listOrgLib, listIntraPic, listRefLib, frameRate, 0);

				% RDcost.
				lambda = lambdaFactor * 2 ^ ((listQP(idxQP) - 12) / 3);
				for idxLibPic = 1: numLibPic
					% for all frames.
					sseY = widthSeq * heightSeq * sseAnchor(idxLibPic, 1);
					bits = bitPsnrAnchor(idxLibPic, 1) / frameRate * listNumFrames(idxLibPic) * 1000;
					cellRDcostAnchorAllLib(idxLibPic, idxQP) = {sseY + lambda * bits};
					% for key frames.
					sseY = widthSeq * heightSeq * sseAnchorKeyFrames(idxLibPic, 1);
					bits = bitPsnrAnchorKeyFrames(idxLibPic, 1) / frameRate * listLibFreq(idxLibPic) * 1000;
					cellRDcostAnchorAllLibKeyFrames(idxLibPic, idxQP) = {sseY + lambda * bits};
				end
			end

			% read libvc.
			% collect data.
			cellRDcostLibvcAllLib = cell(numLibPic, numQP);
			cellRDcostLibvcAllLibKeyFrames = cell(numLibPic, numQP);

			% for each qp.
			for idxQP = 1: numQP
				% prepare log path.
				strQP = num2str(listQP(idxQP));

				% for each dqp.
				% collect data.
				matRDcostALLDqp = zeros(numDqp, numLibPic);
				matRDcostALLDqpKeyFrames = zeros(numDqp, numLibPic);
			
				for idxDqp = 1: numDqp
					% prepare log path.
					if(listDqp(idxDqp) <= 0)
						strDqp = num2str(-listDqp(idxDqp));
					else
						strDqp = ['0' num2str(listDqp(idxDqp))];
					end
					fileLib = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_libvc_gop16_fixQP_' strDqp '_enc.txt'];

					% read libvc log file.
					[bitPsnrLib, bitPsnrLibKeyFrames, sseLib, sseLibKeyFrames, listNumFrames] = readBitPsnrSingleLib(fileLib, listOrgLib, listIntraPic, listRefLib, frameRate, 1);

					% RDcost.
					lambda = lambdaFactor * 2 ^ ((listQP(idxQP) - 12) / 3);
					for idxLibPic = 1: numLibPic
						% for all frames.
						sseY = widthSeq * heightSeq * sseLib(idxLibPic, 1);
						bits = bitPsnrLib(idxLibPic, 1) / frameRate * listNumFrames(idxLibPic) * 1000;
						matRDcostALLDqp(idxDqp, idxLibPic) = sseY + lambda * bits;
						% for key frames.
						sseY = widthSeq * heightSeq * sseLibKeyFrames(idxLibPic, 1);
						bits = bitPsnrLibKeyFrames(idxLibPic, 1) / frameRate * listLibFreq(idxLibPic) * 1000;
						matRDcostALLDqpKeyFrames(idxDqp, idxLibPic) = sseY + lambda * bits;
					end
				end
				% collect data.
				for idxLibPic = 1: numLibPic
					cellRDcostLibvcAllLib(idxLibPic, idxQP) = {matRDcostALLDqp(:, idxLibPic)};
					cellRDcostLibvcAllLibKeyFrames(idxLibPic, idxQP) = {matRDcostALLDqpKeyFrames(:, idxLibPic)};
				end
			end
			% collect data.
			cellRDcostAllSeq(idxSeq, 1) = {cellRDcostAnchorAllLib};
			cellRDcostAllSeq(idxSeq, 2) = {cellRDcostLibvcAllLib};
			cellRDcostAllSeq(idxSeq, 3) = {listNumFrames};

			cellRDcostAllSeqKeyFrames(idxSeq, 1) = {cellRDcostLibvcAllLibKeyFrames};
			cellRDcostAllSeqKeyFrames(idxSeq, 2) = {cellRDcostLibvcAllLibKeyFrames};
			cellRDcostAllSeqKeyFrames(idxSeq, 3) = {listNumFrames};
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			cellRDcostAllSeq_Train = cellRDcostAllSeq;
			cellRDcostAllSeqKeyFrames_Train = cellRDcostAllSeqKeyFrames;
		elseif flagReadSeqListTest
			cellRDcostAllSeq_Test = cellRDcostAllSeq;
			cellRDcostAllSeqKeyFrames_Test = cellRDcostAllSeqKeyFrames;
		end

		% clear temp data.
		clear listIntraPic listOrgLib listRefLib numLibPic widthSeq heightSeq;
		clear cellRDcostAnchorAllLib cellRDcostLibvcAllLib matRDcostALLDqp;
		clear cellRDcostAnchorAllLibKeyFrames cellRDcostLibvcAllLibKeyFrames matRDcostALLDqpKeyFrames;
		clear idxQP idxDqp idxLibPic;
		clear strQP strDqp;
		clear fileAnchor fileLib;
		clear bitPsnrAnchor bitPsnrLib bitPsnrAnchorKeyFrames bitPsnrLibKeyFrames;
		clear listLibPoc listNumFrames;
		clear lambda psnrY sseY bits;
	end
end

if flagPrepareContentInfo
	subFlagComputeContentInfo = 1;
	subsubFlagComputeMC = 1;

    %% read sequence yuv to compute content information of picutures.
	if subFlagComputeContentInfo
		% display information.
		disp('computing content information from sequences.\n');

		oriContentInfoAllSeq = cell(numAllSeq, 3);
        blockSize = 16;
        searchRange = 64;
        
		for idxSeq = 1: numAllSeq
			listIntraPic = cellLibInfo{idxSeq, 4};
			listOrgLib = cellLibInfo{idxSeq, 5};
			listRefLib = cellLibInfo{idxSeq, 6};
			numIntraPic = length(listIntraPic);

			fileYuv = [pathYuv cellListSeq{idxSeq, 1} '\' cellListSeq{idxSeq, 2} '.yuv'];
            disp(fileYuv);

            % \''collect data.
			stdAllIntraPic = zeros(numIntraPic, 1);
			stdBlockAllIntraPic = zeros(numIntraPic, 1);
			distAllIntraPic = zeros(numIntraPic, 1);
			distMCAllIntraPic = zeros(numIntraPic, 2);
            distIntraAllIntraPic = zeros(numIntraPic, 2);
            percentIntraAllIntraPic = zeros(numIntraPic, 1);
            selfIntraAllIntraPic = zeros(numIntraPic, 2);
            storeMv = cell(numIntraPic, 1);
            storeErr = cell(numIntraPic, 1);

            tic;
			parfor idxIntraPic = 1: numIntraPic
				idxLibPic = listRefLib(idxIntraPic) + 1;
                orgPocLibPic = listOrgLib(idxLibPic);
                orgPocIntraPic = listIntraPic(idxIntraPic);
                [yIntraPic] = readYuv(fileYuv, orgPocIntraPic);
                if orgPocLibPic == orgPocIntraPic
                    yLibPic = yIntraPic;
                else
                    [yLibPic] = readYuv(fileYuv, listOrgLib(idxLibPic));
                end		

				% compute variance of luma.
                stdAllIntraPic(idxIntraPic) = std(yIntraPic(:));
                stdBlockAllIntraPic(idxIntraPic) = stdBlock(yIntraPic, 4);

				% compute distance.
				% method1: frame difference.
				distAllIntraPic(idxIntraPic) = mean(abs(yIntraPic(:) - yLibPic(:)));

                % method2: average block difference based on motion compensation.
                if subsubFlagComputeMC
                    if orgPocLibPic ~= orgPocIntraPic
                        [errInterAver, errInterTotal, errIntraAver, errIntraTotal, percentIntra, mvCell, errMatInter] = motionCompensation( yLibPic, yIntraPic, blockSize, searchRange );
                        distMCAllIntraPic(idxIntraPic, :) = [errInterAver errInterTotal];
                        distIntraAllIntraPic(idxIntraPic, :) = [errIntraAver errIntraTotal];
                        percentIntraAllIntraPic(idxIntraPic) = percentIntra;
                        storeMv(idxIntraPic) = {mvCell};
                        storeErr(idxIntraPic) = {errMatInter};
                    else                    
                        [ errIntraAver, errIntraTotal, errMatIntra ] = intraDifference( yLibPic, blockSize );
                        selfIntraAllIntraPic(idxIntraPic, :) = [errIntraAver errIntraTotal];
                        storeErr(idxIntraPic) = {errMatIntra};
                    end
                end
            end
            toc

			oriContentInfoAllSeq(idxSeq, 1) = { [stdAllIntraPic stdBlockAllIntraPic distAllIntraPic distMCAllIntraPic distIntraAllIntraPic percentIntraAllIntraPic selfIntraAllIntraPic] };
			if subsubFlagComputeMC
            	oriContentInfoAllSeq(idxSeq, 2) = {storeMv};
            	oriContentInfoAllSeq(idxSeq, 3) = {storeErr};
            end
		end

		% collect data for train or test.
		if flagReadSeqListTrain
			oriContentInfoAllSeq_Train = oriContentInfoAllSeq;
		elseif flagReadSeqListTest
			oriContentInfoAllSeq_Test = oriContentInfoAllSeq;
		end

		% clear temp data.
		clear listLibFreq listIntraPic listOrgLib listRefLib numIntraPic numLibPic fileYuv;
		clear stdAllIntraPic stdBlockAllIntraPic distAllIntraPic distMCAllIntraPic storeMv storeErr;
		clear yIntraPic uIntraPic vIntraPic idxLibPic yLibPic uLibPic vLibPic;
		clear distMC mvCell errMat;
    end
end

if flagFitForBestDqpEachLib
	subFlagPrepareSampleBestDqpForRDcostEachLib = 1;
	subFlagProcessGainForBestDQPEachLib = 0;
	subFlagProcessContentInfoEachLib = 0;
    subFlagPrepareFitData = 0;

    %% prepare the sample based on RDcost for each library picture.
    if subFlagPrepareSampleBestDqpForRDcostEachLib
    	subsubFlagPrepareSample = 0;
    	subsubFlagTrainFunction = 1;

    	% display information.
		disp('fitting samples.\n');

		if subsubFlagPrepareSample
	    	% collect data.
	    	matSampleAllSeqAllLibAllQp = [];
	    	sample = zeros(1, 6); % SSEi/similarity/numIntra/qstep_key/Dqp_allFrames/qstep_allFrames/Dqp_keyFrames/qstep_keyFrames.

	    	for idxSeq = 1: numAllSeq
	    		% prepare library information.
	    		listLibFreq = cellLibInfo{idxSeq, 3};
	    		listIntraPic = cellLibInfo{idxSeq, 4};
	    		listRefLib = cellLibInfo{idxSeq, 6};
	    		listOrgLib = cellLibInfo{idxSeq, 5};
	    		numLibPic = length(listOrgLib);
	    		numIntraPic = length(listIntraPic);
	    		% prepare content information.
	    		listSimPercentageAllIntraPic = oriContentInfoAllSeq{idxSeq}(:, 8);
	    		listSelfIntraAllIntraPic = oriContentInfoAllSeq{idxSeq}(:, 10); % use the frame-level SSE instead of pixel-level MSE.
	    		widthSeq = cellListSeq{idxSeq, 5};
				heightSeq = cellListSeq{idxSeq, 6};

	    		for idxLibPic = 1: numLibPic
	    			% prepare content information for samples.
	    			posLibPicInIntraList = find(listIntraPic == listOrgLib(idxLibPic));
	    			if isempty(posLibPicInIntraList)
	    				error('error in find the library picture in key picture list\n');
	    			end
	    			sample(1) = listSelfIntraAllIntraPic(posLibPicInIntraList) / widthSeq / heightSeq;

	    			% prepare similarity for each intra picture, find the corresponding library picture and create similarity.
	    			% collect data.
	    			sumSimPercentage = 0;
	    			for idxIntraPic = 1: numIntraPic
	    				% find the corresponding library picture.
	    				if listRefLib(idxIntraPic) ~= (idxLibPic - 1)
	    					continue;
	    				end
	    				% summary the similarity.
	    				sumSimPercentage = sumSimPercentage + listSimPercentageAllIntraPic(idxIntraPic);
	    			end
	    			% since the similarity for the key picture who is library picture is not computed, add the sumsim by 1.
	    			sample(2) = sumSimPercentage + 1;

	    			% prepare number of key pictures for each library picture.
	    			% NOTE!!!!!!!!!!!!!!!--this number should be the total number of pictures that are considered in the RDcost function.
	    			% !!!!!!!!!!!!!!!!!!!--But here we use only the number of key frames.
	    			sample(3) = listLibFreq(idxLibPic);

	    			% prepare QP and best DQP.
	    			for idxQP = 1: numQP
	    				% prepare data.
	    				matRDcostALLDqp = cellRDcostAllSeq{idxSeq, 2}{idxLibPic, idxQP};
	    				matRDcostALLDqpKeyFrames = cellRDcostAllSeqKeyFrames{idxSeq, 2}{idxLibPic, idxQP};

	    				% prepare base q-step for sample.
	    				sample(4) = 2 ^ ((listQP(idxQP) - 4) / 6);

	    				% prepare best DQP for sample.
	    				% all frames.
	    				[~, idxBestDqp] = min(matRDcostALLDqp);
	    				sample(6) = 2 ^ (listDqp(idxBestDqp) / 6); % use q-step instead of DQP.
	                    sample(5) = listDqp(idxBestDqp);
	                    % key frames.
	                    [~, idxBestDqpKeyFrames] = min(matRDcostALLDqpKeyFrames);
	    				sample(8) = 2 ^ (listDqp(idxBestDqpKeyFrames) / 6); % use q-step instead of DQP.
	                    sample(7) = listDqp(idxBestDqpKeyFrames);

	    				% collect sample.
	    				matSampleAllSeqAllLibAllQp = [matSampleAllSeqAllLibAllQp; sample];
	    			end
	    		end
	    	end

	    	% collect data for train or test.
	    	if flagReadSeqListTrain
	    		matSampleAllSeqAllLibAllQp_Train = matSampleAllSeqAllLibAllQp;
	    	elseif flagReadSeqListTest
	    		matSampleAllSeqAllLibAllQp_Test = matSampleAllSeqAllLibAllQp;
	    	end

	    	clear sample;
	    	clear idxSeq idxLibPic idxBestDqp idxBestDqpKeyFrames;
	    	clear listIntraPic listRefLib listOrgLib numLibPic numIntraPic;
	    	clear listSelfIntraAllIntraPic listSimPercentageAllIntraPic;
	    	clear posLibPicInIntraList sumSimPercentage matRDcostALLDqp matRDcostALLDqpKeyFrames;
	    end

    	%% train the derivation of RDcost function based on the sample.
    	if subsubFlagTrainFunction
    		subsubsubFlagLinearRegress = 1;
    		subsubsubFlagNonLinearRegress = 0;

    		% combine multipal input and fit data using linear regress.
    		if subsubsubFlagLinearRegress
                subsubsubsubFlagTestFunction = 1;
                
    			% prepare basic data. SSEi/similarity/numIntra/qstep_key/Dqp_allFrames/qstep_allFrames/Dqp_keyFrames/qstep_keyFrames.
    			sseY = matSampleAllSeqAllLibAllQp_Train(:, 1);
    			sim = matSampleAllSeqAllLibAllQp_Train(:, 2);
    			nIntra = matSampleAllSeqAllLibAllQp_Train(:, 3);
    			qKey = matSampleAllSeqAllLibAllQp_Train(:, 4);
    			%dq = matSampleAllSeqAllLibAllQp_Train(:, 8);
                dq = bestDqpAllSeqAllLibPicAllQp_BDrate_Train;

    			% prepare input data. equation ----- a1*qk^3*n*s*dq^4+a2*qk^2*n*s*dq^4+a3*qk^3*dq^3+a4*qk*sse*dq=sse;
    			% first item.
    			inData_1 = qKey.^3.*nIntra.*sim.*dq.^4;
    			inData_2 = qKey.^2.*nIntra.*sim.*dq.^4;
    			inData_3 = qKey.^3.*dq.^3;
    			inData_4 = qKey.*sseY.*dq;
                inData_5 = ones(length(inData_1), 1);
    			outData = sseY;

    			% linear regress.
    			[aCoeff, ~, ~, ~, stats_Train] = regress(outData, [inData_1 inData_2 inData_3 inData_4 inData_5]);
    			aCoeff
    			stats_Train(1)

    			if subsubsubsubFlagTestFunction
    				
    				% prepare input data.
                    matSample = matSampleAllSeqAllLibAllQp_Train;
    				dq_FitBest = zeros(size(matSample, 1), 1) - 100;
                    valMin = Inf(size(matSample, 1), 1);
    				for idxSample = 1: length(dq_FitBest)
    					% prepare basic data.
	    				sseY_Test = matSample(idxSample, 1);
	    				sim_Test = matSample(idxSample, 2);
	    				nIntra_Test = matSample(idxSample, 3);
	    				qKey_Test = matSample(idxSample, 4);
	    				% check for the best.
    					for idxDqp = 1: numDqp
    						dq_FitTmp = 2 ^ (listDqp(idxDqp) / 6);
    						% prepare data.
    						inData_1 = qKey_Test.^3.*nIntra_Test.*sim_Test.*dq_FitTmp.^4;
			    			inData_2 = qKey_Test.^2.*nIntra_Test.*sim_Test.*dq_FitTmp.^4;
			    			inData_3 = qKey_Test.^3.*dq_FitTmp.^3;
			    			inData_4 = qKey_Test.*sseY_Test.*dq_FitTmp;
                            inData_5 = ones(length(inData_1), 1);
			    			outData = sseY_Test;

    						valTmp = aCoeff(1)*inData_1 ...
    								+ aCoeff(2)*inData_2 ...
    								+ aCoeff(3)*inData_3 ...
    								+ aCoeff(4)*inData_4 ...
                                    + aCoeff(5)*inData_5 ...
    								- outData;
    						if abs(valTmp) < valMin(idxSample)
    							dq_FitBest(idxSample) = listDqp(idxDqp);
    							valMin(idxSample) = abs(valTmp);
    						end
    					end
    				end
    				% clear temp data.
    				clear idxSample idxDqp;
    				clear valTmp;
    				clear sseY_Test sim_Test nIntra_Test qKey_Test;
    				clear dq_FitTmp;
    			end

    			clear sseY sim nIntra qKey dq;
    			clear inData_1 inData_2 inData_3 inData_4 inData_5 outData;
    		end

			% combine multipal input and fit data using linear regress.
    		if subsubsubFlagNonLinearRegress
				% prepare basic data. SSEi/similarity/numIntra/qstep_key/Dqp_allFrames/qstep_allFrames/Dqp_keyFrames/qstep_keyFrames.
    			sseY = matSampleAllSeqAllLibAllQp_Train(:, 1);
    			sim = matSampleAllSeqAllLibAllQp_Train(:, 2);
    			nIntra = matSampleAllSeqAllLibAllQp_Train(:, 3);
    			qKey = matSampleAllSeqAllLibAllQp_Train(:, 4);
    			dq = matSampleAllSeqAllLibAllQp_Train(:, 8);
    			outData = zeros(length(dq), 1);

				% Set up fittype and options. SSEi/similarity/numIntra/qstep_key/Dqp_allFrames/qstep_allFrames/Dqp_keyFrames/qstep_keyFrames.
				equation = @(a,x) ((a(1)*x(4)^3+a(2)*x(4)^2)*x(2)*x(3)*x(5)^4+a(3)*x(4)^3*x(5)^3-a(4)*x(4)*x(1)*x(5)-a(5)*x(5));
				
				% Fit model to data.
				aInit = [1;1;1;1;1];
				[aCoeff,R,J,CovB,MSE,ErrorModelInfo] = nlinfit([sseY sim nIntra qKey dq], outData, equation, aInit);
			end
    	end
    end

	%% compute the best non-integer delta QP with curve fit.
	if subFlagProcessGainForBestDQPEachLib
		% fetch the gain data, which can be modified.
		gainAllSeq = cell2mat(gainYEachLibAllSeq);

		% collect data.
		bestNonIntDqpAllSeq = zeros(size(gainAllSeq, 1), 1);

		% prepare for fit.
		xDqp = listDqp;
		%equation = 'a*(x-b)^2+c';
		equation = 'a*x+b/(x-3)+c';
        countTemp = 0;

		% operate curve fit for all sequence to obtain best DQP.
		for idxSeq = 1: numAllSeq
			% for each library picture.
			numLibPic = length(cellLibInfo{idxSeq, 5});
			for idxLibPic = 1: numLibPic
				countTemp = countTemp + 1;
				gain = gainAllSeq(countTemp, :);

				% fit curve.
                startPoint = [1 1 1];
                if countTemp == 11
                    startPoint = [0.001 1 0.001];
                end
				[ fitresult, gof ] = computeEquationFit( xDqp, gain, equation, startPoint, 0 );
                gof
				% best non-integer DQP.
				bestNonIntDqpAllSeq(countTemp) = fitresult.b;
			end
		end

		% clear temp data.
		clear gainAllSeq;
		clear xDqp equation countTemp;
		clear numLibPic gain fitresult gof;
	end

	%% process the content information for model fitting.
	if subFlagProcessContentInfoEachLib
		% collect data for all sequences.
		numEachLib = size(bestNonIntDqpAllSeq, 1);
		contentInfoEachLib = zeros(numEachLib, 3, 3);

		% prepare for loop.
		countTemp = 0;

		for idxSeq = 1: numAllSeq
			% prepare library information.
			listLibFreq = cellLibInfo{idxSeq, 3}';
			listOrgLib = cellLibInfo{idxSeq, 5}';
			listRefLib = cellLibInfo{idxSeq, 6}';
			numLibPic = length(listOrgLib);

			% prepare content information.
			stdAllIntraPic = oriContentInfoAllSeq{idxSeq, 1}(:, 1);
			stdBlockAllIntraPic = oriContentInfoAllSeq{idxSeq, 1}(:, 2);
            distAllIntraPic = oriContentInfoAllSeq{idxSeq, 1}(:, 3);
            distMCAllIntraPic = oriContentInfoAllSeq{idxSeq, 1}(:, 4);

			% collect content information for each library.
			averTextStd = zeros(numLibPic, 1);
			averTextStdBlock = zeros(numLibPic, 1);
			sumDist = zeros(numLibPic, 1);
			sumDistMC = zeros(numLibPic, 1);
            sumDistLog = zeros(numLibPic, 1);

			% process for each library picture.
			for idxLibPic = 1: numLibPic
				% prepare to get all key frame referencing the library picture.
				idxKeyPicRefLibPic = find(listRefLib == idxLibPic - 1);
				numKeyPicRefLibPic = length(idxKeyPicRefLibPic);
				if numKeyPicRefLibPic ~= listLibFreq(idxLibPic)
					error('number of key pictures referencing one library picture is not equal to the frequency.\n');
				end

				% compute summary texture.
				% method 1: texture based on std of whole picture.
				averTextStd(idxLibPic) = sum(stdAllIntraPic(idxKeyPicRefLibPic)) / numKeyPicRefLibPic;
				% method 2: texture based on average std of all blocks in a picture.
				averTextStdBlock(idxLibPic) = sum(stdBlockAllIntraPic(idxKeyPicRefLibPic)) / numKeyPicRefLibPic;

				% compute summary distance.
				% method 1: frame difference.
				sumDist(idxLibPic) = sum(distAllIntraPic(idxKeyPicRefLibPic));
				% method 2: difference based on motion compensation.
				sumDistMC(idxLibPic) = sum(distMCAllIntraPic(idxKeyPicRefLibPic));
                % method 3: log frame difference.
                for i = 1: length(idxKeyPicRefLibPic)
                    if distAllIntraPic(idxKeyPicRefLibPic(i)) == 0
                        continue;
                    end
                    sumDistLog(idxLibPic) = sumDistLog(idxLibPic) + log(distAllIntraPic(idxKeyPicRefLibPic(i)));
                end

				% collect data.
				countTemp = countTemp + 1;
				contentInfoEachLib(countTemp, 1, 1) = averTextStd(idxLibPic);          
				contentInfoEachLib(countTemp, 2, 1) = averTextStdBlock(idxLibPic);
				contentInfoEachLib(countTemp, 1, 2) = sumDist(idxLibPic);
	            contentInfoEachLib(countTemp, 2, 2) = sumDistMC(idxLibPic);
                contentInfoEachLib(countTemp, 3, 2) = sumDistLog(idxLibPic);
	            contentInfoEachLib(countTemp, 1, 3) = numKeyPicRefLibPic;
			end
		end

		% clear temp data.
		clear numEachLib countTemp;
		clear listLibFreq listIntraPic listOrgLib listRefLib numIntraPic numLibPic intraPeriod;
		clear stdAllIntraPic stdBlockAllIntraPic distAllIntraPic distMCAllIntraPic mvAllIntraPic;
		clear averTextStd sumDist sumDistMC sumDistLog;
		clear idxKeyPicRefLibPic numKeyPicRefLibPic;
    end
    
    %% prepare data for curve fit.
    if subFlagPrepareFitData
        
    end
end

%% check the result.
if flagCheckResultValidity
	subFlagFindBestDqpForAllSeqAllLibPicAllQp = 1;

	%% find the best delta QP for all sequences, all library pictures, all QP.
    if subFlagFindBestDqpForAllSeqAllLibPicAllQp
    	% display information.
		disp('finding the best delta QP for all sequences, all library pictures, all QP.\n');
		% collect data.
		bestDqpAllSeqAllLibPicAllQp_BDrate = [];
		bestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate = [];
        bestGainAllSeqAllLibPic_BDrate = [];
        bestGainAllSeqAllLibPicKeyFrames_BDrate = [];

		for idxSeq = 1: numAllSeq
	    	% prepare data.
			listIntraPic = cellLibInfo{idxSeq, 4};
			listOrgLib = cellLibInfo{idxSeq, 5};
			listRefLib = cellLibInfo{idxSeq, 6};
	    	numLibPic = length(listOrgLib);
	        frameRate = cellListSeq{idxSeq, 3};
			% read anchor.
			% collect data.
			bitPsnrAnchorAllLib = cell(numQP, numLibPic);
			bitPsnrAnchorAllLibKeyFrames = cell(numQP, numLibPic);
			for idxQP = 1: numQP
				strQP = num2str(listQP(idxQP));
				fileAnchor = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_anchor_gop16_enc.txt'];
				[bitPsnrAnchor, bitPsnrAnchorKeyFrames, ~, ~, listNumFrames] = readBitPsnrSingleLib(fileAnchor, listOrgLib, listIntraPic, listRefLib, frameRate, 0);
				for idxLibPic = 1: numLibPic
					bitPsnrAnchorAllLib(idxQP, idxLibPic) = {bitPsnrAnchor(idxLibPic, 1: 4)};
					bitPsnrAnchorAllLibKeyFrames(idxQP, idxLibPic) = {bitPsnrAnchorKeyFrames(idxLibPic, 1: 4)};
				end
			end

			% collect data.
			bitPsnrLibPicAllLibAllDqp = cell(numQP, numDqp, numLibPic);
			bitPsnrLibPicAllLibAllDqpKeyFrames = cell(numQP, numDqp, numLibPic);
			for idxDqp = 1: numDqp
				% prepare path.
				if(listDqp(idxDqp) <= 0)
					strDqp = num2str(-listDqp(idxDqp));
				else
					strDqp = ['0' num2str(listDqp(idxDqp))];
				end
				% collect data.
				for idxQP = 1: numQP
					strQP = num2str(listQP(idxQP));
					fileLib = [path cellListSeq{idxSeq, 2} '\Sequence_' strQP '_libvc_gop16_fixQP_' strDqp '_enc.txt'];
					[bitPsnrLib, bitPsnrLibKeyFrames, ~, ~, listNumFrames] = readBitPsnrSingleLib(fileLib, listOrgLib, listIntraPic, listRefLib, frameRate, 1);
					for idxLibPic = 1: numLibPic
						bitPsnrLibPicAllLibAllDqp(idxQP, idxDqp, idxLibPic) = {bitPsnrLib(idxLibPic, 1: 4)};
						bitPsnrLibPicAllLibAllDqpKeyFrames(idxQP, idxDqp, idxLibPic) = {bitPsnrLibKeyFrames(idxLibPic, 1: 4)};
					end
				end
			end

			% compute BD-rate and find the dest delta QP.
			for idxLibPic = 1: numLibPic
				% BD-rate for all frames.
				bestGain = Inf;
				bestGainKeyFrames = Inf;
				bestDqpIdx = zeros(4, 1);
				bestDqpIdxKeyFrames = zeros(4, 1);
				% temporal memory for QP 22, 27, 32, 37.
				bitPsnrAnchor = [bitPsnrAnchorAllLib{1, idxLibPic}; bitPsnrAnchorAllLib{2, idxLibPic}; bitPsnrAnchorAllLib{3, idxLibPic}; bitPsnrAnchorAllLib{4, idxLibPic};];
				bitPsnrAnchorKeyFrames = [bitPsnrAnchorAllLibKeyFrames{1, idxLibPic}; bitPsnrAnchorAllLibKeyFrames{2, idxLibPic}; bitPsnrAnchorAllLibKeyFrames{3, idxLibPic}; bitPsnrAnchorAllLibKeyFrames{4, idxLibPic};];
				bitPsnrLib = zeros(4, 4);
				bitPsnrLibKeyFrames = zeros(4, 4);

				for idxDqp1 = 1: numDqp
					bitPsnrLib(1, :) = bitPsnrLibPicAllLibAllDqp{1, idxDqp1, idxLibPic};
					bitPsnrLibKeyFrames(1, :) = bitPsnrLibPicAllLibAllDqpKeyFrames{1, idxDqp1, idxLibPic};
					for idxDqp2 = 1: numDqp
						bitPsnrLib(2, :) = bitPsnrLibPicAllLibAllDqp{2, idxDqp2, idxLibPic};
						bitPsnrLibKeyFrames(2, :) = bitPsnrLibPicAllLibAllDqpKeyFrames{2, idxDqp2, idxLibPic};
						for idxDqp3 = 1: numDqp
							bitPsnrLib(3, :) = bitPsnrLibPicAllLibAllDqp{3, idxDqp3, idxLibPic};
							bitPsnrLibKeyFrames(3, :) = bitPsnrLibPicAllLibAllDqpKeyFrames{3, idxDqp3, idxLibPic};
							for idxDqp4 = 1: numDqp
								bitPsnrLib(4, :) = bitPsnrLibPicAllLibAllDqp{4, idxDqp4, idxLibPic};
								bitPsnrLibKeyFrames(4, :) = bitPsnrLibPicAllLibAllDqpKeyFrames{4, idxDqp4, idxLibPic};

								% compute BD-rate.
								gainTmp = bdRateComparation( bitPsnrAnchor, bitPsnrLib );
								if gainTmp(1) < bestGain
									bestDqpIdx = [idxDqp1; idxDqp2; idxDqp3; idxDqp4];
                                    bestGain = gainTmp(1);
								end
								gainTmpKeyFrames = bdRateComparation( bitPsnrAnchorKeyFrames, bitPsnrLibKeyFrames );
								if gainTmpKeyFrames < bestGainKeyFrames
									bestDqpIdxKeyFrames = [idxDqp1; idxDqp2; idxDqp3; idxDqp4];
                                    bestGainKeyFrames = gainTmpKeyFrames(1);
								end
							end
						end
					end
				end
				% collect data.
				bestDqpAllSeqAllLibPicAllQp_BDrate = [bestDqpAllSeqAllLibPicAllQp_BDrate; listDqp(bestDqpIdx)'];
				bestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate = [bestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate; listDqp(bestDqpIdxKeyFrames)'];
				bestGainAllSeqAllLibPic_BDrate = [bestGainAllSeqAllLibPic_BDrate; bestGain];
				bestGainAllSeqAllLibPicKeyFrames_BDrate = [bestGainAllSeqAllLibPicKeyFrames_BDrate; bestGainKeyFrames];
			end
		end

		if flagReadSeqListTrain
			bestDqpAllSeqAllLibPicAllQp_BDrate_Train = bestDqpAllSeqAllLibPicAllQp_BDrate;
			bestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate_Train = bestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate;
			bestGainAllSeqAllLibPic_BDrate_Train = bestGainAllSeqAllLibPic_BDrate;
			bestGainAllSeqAllLibPicKeyFrames_BDrate_Train = bestGainAllSeqAllLibPicKeyFrames_BDrate;
		elseif flagReadSeqListTest
			bestDqpAllSeqAllLibPicAllQp_BDrate_Test = bestDqpAllSeqAllLibPicAllQp_BDrate;
			bestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate_Test = bestDqpAllSeqAllLibPicAllQpKeyFrames_BDrate;
			bestGainAllSeqAllLibPic_BDrate_Test = bestGainAllSeqAllLibPic_BDrate;
			bestGainAllSeqAllLibPicKeyFrames_BDrate_Test = bestGainAllSeqAllLibPicKeyFrames_BDrate;
		end

		% clear temp data.
		clear frameRate dqp qp strDqp strQP fileLibInfo numLibPic strDqp fileLib;
		clear strQP fileAnchor listLibPoc listNumFrames;
		clear gain bitPsnrAnchor bitPsnrLib bitPsnrAnchorKeyFrames bitPsnrLibKeyFrames;
		clear bitPsnrAnchorAllLib bitPsnrLibAllLib bitPsnrAnchorAllLibKeyFrames bitPsnrLibAllLibKeyFrames;
		clear bestGain bestGainKeyFrames bestDqpIdx bestDqpIdxKeyFrames gainTmp gainTmpKeyFrames ;
		clear idxDqp1 idxDqp2 idxDqp3 idxDqp4 idxDqp idxLibPic;
		clear bitPsnrLibPicAllLibAllDqp bitPsnrLibPicAllLibAllDqpKeyFrames;
    end
end