%% This script is used to read library information for all sequences with all QPs and delta QPs.
flagReadSeqListCTC = 0;
flagReadSeqListTrain = 0;
flagReadLibInfo = 0;
flagReadLogFile = 0;
flagPrepareContentInfo = 1;
flagFitForBestDqpEachLib = 0;
flagFitForBestDqpWholeSeq = 0;

listQP = 22:5:37;
listDqp = 2:-1:-11;
numQP = length(listQP);
numDqp = length(listDqp);

%% read sequences list file of CTC.
if flagReadSeqListCTC
    rcVecSeqsB = readSeqList('.\\SeqList_B.txt');
    rcVecSeqsC = readSeqList('.\\SeqList_C.txt');
    rcVecSeqsD = readSeqList('.\\SeqList_D.txt');
    numVecSeqsB = size(rcVecSeqsB, 1);
    numVecSeqsC = size(rcVecSeqsC, 1);
    numVecSeqsD = size(rcVecSeqsD, 1);
    classSeqB = cellstr(char(ones(numVecSeqsB, 1) * 'B'));
    classSeqC = cellstr(char(ones(numVecSeqsC, 1) * 'C'));
    classSeqD = cellstr(char(ones(numVecSeqsD, 1) * 'D'));

    % collect data.
    cellListSeq = [classSeqB rcVecSeqsB; classSeqC rcVecSeqsC ; classSeqD rcVecSeqsD];
    numAllSeq = numVecSeqsB + numVecSeqsC + numVecSeqsD;

    % set path.
    path = '.\log\';
    pathYuv = '\\mcl\MCL_Space\TestSeq\FinalCfPVersionSequence\';
    
    % clear temp data.
    clear rcVecSeqsB rcVecSeqsC rcVecSeqsD numVecSeqsB numVecSeqsC numVecSeqsD classSeqB classSeqC classSeqD;
end
%% read sequences list file of AVS.
if flagReadSeqListTrain
    rcVecSeqs = readSeqList('.\\SeqList-train-all.txt');
    numVecSeqs = size(rcVecSeqs, 1);

    % collect data.
    cellListSeq = rcVecSeqs;
    numAllSeq = numVecSeqs;
    
    % set path.
    path = '.\log-train\';
    pathYuv = '\\mcl\MCL_Space\TestSeq\';

    % clear temp data.
    clear rcVecSeqs numVecSeqs;
end

%% read library information file.
if flagReadLibInfo
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

	% clear temp data.
	clear dqp qp strDqp strQP fileLibInfo libInfo listSCD listOrgLib listIntraPic listRefLib listLibFreq listCenterCost listClusterCost;
end

if flagReadLogFile
	subFlagReadSeqLog = 0;
	subFlagReadKeyLog = 0;
	subFlagReadBitPsnrEachLib = 1;

	%% read sequence encoding log file for bitrate and psnr.
	if subFlagReadSeqLog
		gainYUVAllSeq = zeros(numAllSeq, numDqp);
	    bitPsnrALLSeqAnchor = cell(numAllSeq, 1);
		for idxSeq = 1: numAllSeq
			% read anchor.
			bitPsnrAnchor = zeros(4, 5);
			for idxQP = 1: numQP
				strQP = num2str(listQP(idxQP));
				fileAnchor = [path cellListSeq{idxSeq, 2} '\seq\Sequence_' strQP '_anchor_enc.txt'];
				bitPsnrAnchor(idxQP, :) = readSumBitPsnr(fileAnchor);
			end
			bitPsnrALLSeqAnchor(idxSeq) = {bitPsnrAnchor};

			for idxDqp = 1: numDqp
				if(listDqp(idxDqp) <= 0)
					strDqp = num2str(-listDqp(idxDqp));
				else
					strDqp = ['0' num2str(listDqp(idxDqp))];
				end
				bitPsnrLib = zeros(4, 5);
				for idxQP = 1: numQP
					strQP = num2str(listQP(idxQP));
					fileLib = [path cellListSeq{idxSeq, 2} '\seq\Sequence_fixQP_' strDqp '_' strQP '_libvc_enc.txt'];
					bitPsnrLib(idxQP, :) = readSumBitPsnr(fileLib);
				end
				gain = bdRateComparation( bitPsnrAnchor, bitPsnrLib );
				gainYUVAllSeq(idxSeq, idxDqp) = gain(4);
			end
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
				fileAnchor = [path cellListSeq{idxSeq, 2} '\seq\Sequence_' strQP '_anchor_enc.txt'];
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
					fileLib = [path cellListSeq{idxSeq, 2} '\seq\Sequence_fixQP_' strDqp '_' strQP '_libvc_enc.txt'];
					bitPsnrLib(idxQP, :) = readKeyFrameSumBitPsnr(fileLib);
				end
				gain = bdRateComparation( bitPsnrAnchor, bitPsnrLib );
				gainKeyYUVAllSeq(idxSeq, idxDqp) = gain(4);
			end
		end
		% clear temp data.
		clear bitPsnrAnchor strQP fileAnchor strDqp bitPsnrLib fileLib gain;
	end

	%% read average bitrate and psnr of frames for each library picture.
	if subFlagReadBitPsnrEachLib
		gainYEachLibAllSeq = cell(numAllSeq, numDqp);
	    bitPsnrEachLibALLSeqLibvc = cell(numAllSeq, numDqp);
	    numFrameEachLibAllSeq = cell(numAllSeq, 1);
		for idxSeq = 1: numAllSeq
	        frameRate = cellListSeq{idxSeq, 4} / 10;
			% for one sequence, libInfo is same for all qp and dqp, except the area of inter prediction.
	    	dqp = 0; qp = 22;
	    	strDqp = num2str(dqp);
	    	strQP = num2str(qp);
	    	fileLibInfo = [path cellListSeq{idxSeq, 2} '\lib\libraryInfo_fixQP_' strDqp '_' strQP '_libvc.txt'];
	    	numLibPic = length(cellLibInfo{idxSeq, 5});
			% read anchor.
			bitPsnrAnchorAllLib = zeros(4* numLibPic, 5);
			for idxQP = 1: numQP
				strQP = num2str(listQP(idxQP));
				fileAnchor = [path cellListSeq{idxSeq, 2} '\seq\Sequence_' strQP '_anchor_enc.txt'];
				[bitPsnrAnchor, listLibPoc, listNumFrames] = readBitPsnrSingleLib(fileAnchor, fileLibInfo, frameRate, 0);
				bitPsnrAnchorAllLib(idxQP: 4: end, :) = bitPsnrAnchor;
	            numFrameEachLibAllSeq(idxSeq) = {listNumFrames};
			end
			numFrameEachLibAllSeq(idxSeq) = {listNumFrames};

			for idxDqp = 1: numDqp
				if(listDqp(idxDqp) <= 0)
					strDqp = num2str(-listDqp(idxDqp));
				else
					strDqp = ['0' num2str(listDqp(idxDqp))];
				end
				bitPsnrLibAllLib = zeros(4* numLibPic, 5);
				for idxQP = 1: numQP
					strQP = num2str(listQP(idxQP));
					fileLib = [path cellListSeq{idxSeq, 2} '\seq\Sequence_fixQP_' strDqp '_' strQP '_libvc_enc.txt'];
					[bitPsnrLib, listLibPoc, listNumFrames] = readBitPsnrSingleLib(fileLib, fileLibInfo, frameRate, 1);
					bitPsnrLibAllLib(idxQP: 4: end, :) = bitPsnrLib;
				end
				bitPsnrEachLibALLSeqLibvc(idxSeq, idxDqp) = {bitPsnrLibAllLib};

				% compute BD-rate for each group of library picture.
				gain = zeros(numLibPic, 4);
				for idxLibPic = 1: numLibPic
					bitPsnrAnchor = bitPsnrAnchorAllLib(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
					bitPsnrLib = bitPsnrLibAllLib(((idxLibPic - 1) * 4 + 1): (idxLibPic * 4), :);
					gain(idxLibPic, :) = bdRateComparation( bitPsnrAnchor, bitPsnrLib );
				end
				gainYEachLibAllSeq(idxSeq, idxDqp) = {gain(:, 1)};
			end
		end
		% clear temp data.
		clear frameRate dqp qp strDqp strQP fileLibInfo numLibPic bitPsnrAnchorAllLib;
		clear strQP fileAnchor listLibPoc listNumFrames;
		clear strDqp bitPsnrLibAllLib fileLib;
		clear gain bitPsnrAnchor bitPsnrLib;
	end
end

if flagPrepareContentInfo
	subFlagComputeContentInfo = 1;
	subsubFlagComputeMC = 1;

    %% read sequence yuv to compute content information of picutures.
	if subFlagComputeContentInfo
		
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
		% clear temp data.
		clear listLibFreq listIntraPic listOrgLib listRefLib numIntraPic numLibPic fileYuv;
		clear stdAllIntraPic stdBlockAllIntraPic distAllIntraPic distMCAllIntraPic storeMv storeErr;
		clear yIntraPic uIntraPic vIntraPic idxLibPic yLibPic uLibPic vLibPic;
		clear distMC mvCell errMat;
    end
end

if flagFitForBestDqpEachLib
	subFlagProcessGainForBestDQPEachLib = 0;
	subFlagProcessContentInfoEachLib = 1;
    subFlagPrepareFitData = 0;

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
        subsubFlagQStep = 1;
        if subsubFlagQStep
            zTr = 2 .^ (bestNonIntDqpAllSeq_Train / 6);
            zTe = 2 .^ (bestNonIntDqpAllSeq_Test / 6);
        else
            zTr = bestNonIntDqpAllSeq_Train;
            zTe = bestNonIntDqpAllSeq_Test;
        end
        zTr = zTr .^ 2;
        zTe = zTe .^ 2;
        
        idxCase = [5 13];
        
        pageStd = 1; pageDist = 2; pageFreq = 3;
        % data index in order of std, dist, freq.
        idxData = [
            1 1 1;    % 1--std, sumDist, freq.
            1 2 1;    % 2--std, sumDistMC, freq.        x
            2 1 1;    % 3--stdBlock, sumDist, freq.
            2 2 1;    % 4--stdBlock, sumDistMC, freq.   x
            1 3 1;    % 5--std, sumDistLog, freq.
            2 3 1;    % 6--stdBlock, sumDistLog, freq.
            ];
        % train data.
        dataStd = contentInfoEachLib_Train(:, idxData(idxCase(1), pageStd), pageStd);
        dataDist = contentInfoEachLib_Train(:, idxData(idxCase(1), pageDist), pageDist);
        dataFreq = contentInfoEachLib_Train(:, idxData(idxCase(1), pageFreq), pageFreq);
        dataTr = [
            % 1--x = (std)/(sumDist);y = (freq)/(std);              x
            {dataStd./dataDist} {dataFreq./dataStd};
            % 2--x = (sumDist)/(std);y = (std)/(freq);              x
            {dataDist./dataStd} {dataStd./dataFreq};
            % 3--x = (averDist)/(std);y = (std)/(freq);       
            {dataDist./dataStd./dataFreq} {dataStd./dataFreq};
            % 4--x = (std)/(averDist);y = (freq)/(std);
            {dataStd./dataDist.*dataFreq} {dataFreq./dataStd};
            % 5--x = 1/(averDist*std);y = (freq)/(std);             x
            {1./dataDist.*dataFreq./dataStd} {dataFreq./dataStd};
            % 6--x = averDist*std;y = (freq)/(std);                 x
            {dataDist./dataFreq.*dataStd} {dataFreq./dataStd};
            % 7--x = averDist/std;y = (1/freq);                       x
            {dataDist./dataFreq.*dataStd} {1./dataFreq};
            % 8--x = averDist;y = (1/freq);
            {dataDist./dataFreq} {1./dataFreq};
            % 9--x = log(averDist);y = (1/log(freq))
            {log(dataDist./dataFreq)} {1./log(dataFreq)};
            % 10--x = sumdist; y=freq
            {dataDist./dataFreq} {dataFreq}
            % 11--x = sumdist / freq; y=std/freq
            {dataDist./dataFreq} {dataStd./dataFreq}
            % 12--x = sumdist * std / freq; y=std/freq
            {dataDist.*dataStd./dataFreq} {dataStd./dataFreq}
            % 13--x = sumdist * std / freq; y=log(std)/freq
            {dataDist.*log2(dataStd)./dataFreq} {log2(dataStd)./dataFreq}
            ];
        % test data.
        dataStd = contentInfoEachLib_Test(:, idxData(idxCase(1), pageStd), pageStd);
        dataDist = contentInfoEachLib_Test(:, idxData(idxCase(1), pageDist), pageDist);
        dataFreq = contentInfoEachLib_Test(:, idxData(idxCase(1), pageFreq), pageFreq);
        dataTe = [
            % 1--x = (std)/(sumDist);y = (freq)/(std);
            {dataStd./dataDist} {dataFreq./dataStd};
            % 2--x = (sumDist)/(std);y = (std)/(freq);
            {dataDist./dataStd} {dataStd./dataFreq};
            % 3--x = (averDist)/(std);y = (std)/(freq);
            {dataDist./dataStd./dataFreq} {dataStd./dataFreq};
            % 4--x = (std)/(averDist);y = (freq)/(std);
            {dataStd./dataDist.*dataFreq} {dataFreq./dataStd};
            % 5--x = 1/(averDist*std);y = (freq)/(std);
            {1./dataDist.*dataFreq./dataStd} {dataFreq./dataStd};
            % 6--x = averDist*std;y = (freq)/(std);
            {dataDist./dataFreq.*dataStd} {dataFreq./dataStd};
           % 7--x = averDist/std;y = (1/freq);                       x
            {dataDist./dataFreq.*dataStd} {1./dataFreq};
            % 8--x = averDist;y = (1/freq);
            {dataDist./dataFreq} {1./dataFreq};
            % 9--x = log(averDist);y = (1/log(freq))
            {log(dataDist./dataFreq)} {1./log(dataFreq)};
            % 10--x = sumdist; y=freq
            {dataDist./dataFreq} {dataFreq}
            % 11--x = sumdist / freq; y=std/freq
            {dataDist./dataFreq} {dataStd./dataFreq}
            % 12--x = sumdist * std / freq; y=std/freq
            {dataDist.*dataStd./dataFreq} {dataStd./dataFreq}
            % 13--x = sumdist * log(std) / freq; y=log(std)/freq
            {dataDist.*log2(dataStd)./dataFreq} {log2(dataStd)./dataFreq}
            ];

        xTr = dataTr{idxCase(2), 1};
        yTr = dataTr{idxCase(2), 2};
        xTe = dataTe{idxCase(2), 1};
        yTe = dataTe{idxCase(2), 2};
        
        [ fitresult, gof, zTeSim, rmseTe ] = computeEquationSurfaceFit( xTr, yTr, zTr, xTe, yTe, zTe, 'poly22', 1 );
        %[ fitresult, gof, zTeSim, rmseTe ] = computeEquationSurfaceFit( xTe, yTe, zTe, xTr, yTr, zTr, 'poly22', 1 );
        
        zTeSim = zTeSim.^0.5;
        if subsubFlagQStep
            zTeSim(zTeSim<0.2) = 0.2;
            zTeSim = log2(zTeSim) * 6;
        end
        zTeSimInt = round(zTeSim);
        zTeSimInt(zTeSimInt<-11) = -11;
        zTeSimInt(zTeSimInt>2) = 2;

        % compute corresponding BD-rate.
        gainTotal = cell2mat(gainYEachLibAllSeq_Test);
        [~, listBestDqp] = min(gainTotal, [], 2);
        listBestDqp = 3 - listBestDqp;
        listSimDqp = zTeSimInt;
        gainMax = computeTotalGain(bitPsnrEachLibALLSeqLibvc_Test, bitPsnrALLSeqAnchor_Test, numFrameEachLibAllSeq_Test, listBestDqp);
        gainSim = computeTotalGain(bitPsnrEachLibALLSeqLibvc_Test, bitPsnrALLSeqAnchor_Test, numFrameEachLibAllSeq_Test, listSimDqp);
        
        output = [listSimDqp;gof.rsquare; gof.rmse; rmseTe; gainSim(:, 1)];
    end
end

%% process library information.
if flagFitForBestDqpWholeSeq
	subFlagProcessContentInfo = 1;
	subFlagBestQP = 0;
    subFlagBestQPKeyFrame = 1;
    subFlagShowData = 1;

    %% process the content information for model fitting.
	if subFlagProcessContentInfo
		
		contentInfoAllSeq = zeros(numAllSeq, 7, 2);

		for idxSeq = 1: numAllSeq
			% mutual information.
	    	listClusterCost = cellLibInfo{idxSeq, 2};
            listClusterCost = listClusterCost';

			% simple content information.
			listLibFreq = cellLibInfo{idxSeq, 3};
			listIntraPic = cellLibInfo{idxSeq, 4};
			listOrgLib = cellLibInfo{idxSeq, 5};
			listRefLib = cellLibInfo{idxSeq, 6};
			numIntraPic = length(listIntraPic);
			numLibPic = length(listOrgLib);
			intraPeriod = cellListSeq{idxSeq, 3};

            averTextureStd =0;
            stdAllIntraPic = oriContentInfoAllSeq{idxSeq, 1}(:, 2);
            distAllIntraPic = oriContentInfoAllSeq{idxSeq, 1}(:, 3);
            distMCAllIntraPic = oriContentInfoAllSeq{idxSeq, 1}(:, 4);
            mvAllIntraPic = oriContentInfoAllSeq{idxSeq, 2}(:, 1);

			sumDistAllLibPic = zeros(numLibPic, 1);
            sumSimAllLibPic = zeros(numLibPic, 1);

            sumDistPocAllLibPic = zeros(numLibPic, 1);
            sumSimPocAllLibPic = zeros(numLibPic, 1);
            
            sumDistMCAllLibPic = zeros(numLibPic, 1);
            sumEnergyMvAllLibPic = zeros(numLibPic, 1);
            
			for idxIntraPic = 1: numIntraPic
				idxLibPic = listRefLib(idxIntraPic) + 1;

				% calculate texture. As for texture, all pictures including key frames that are library pictures should be considered.
				% method 1: variance of luma.
                averTextureStd = averTextureStd + stdAllIntraPic(idxIntraPic);

                % skip if the intra picture is the library picture.
                if listOrgLib(idxLibPic) == listIntraPic(idxIntraPic)
                	continue;
                end
								
				% compute distance.
				% method 1: frame difference.
				dist = distAllIntraPic(idxIntraPic);
				if dist == 0
					err('the distance is zeros\n');
				end

				sumDistAllLibPic(idxLibPic) = sumDistAllLibPic(idxLibPic) + dist;
                sumSimAllLibPic(idxLibPic) = sumSimAllLibPic(idxLibPic) + 1 / dist;

                % method 2: frame poc difference.
                distPoc = abs(listIntraPic(idxIntraPic) - listOrgLib(idxLibPic)) / intraPeriod;
                if distPoc == 0
					err('the distance is zeros\n');
				end
                sumDistPocAllLibPic(idxLibPic) = sumDistPocAllLibPic(idxLibPic) + distPoc;
                sumSimPocAllLibPic(idxLibPic) = sumSimPocAllLibPic(idxLibPic) + 1 / distPoc;
                
                % method 3: difference based on motion compensation.
                distMC = distMCAllIntraPic(idxIntraPic);
                sumDistMCAllLibPic(idxLibPic) = sumDistMCAllLibPic(idxLibPic) + distMC;

                % method 4: energy of mv.
                mvCell = mvAllIntraPic{idxIntraPic};
                [rMvCell, cMvCell] = size(mvCell);
                energyMv = 0;
                for idxR = 1: rMvCell
                	for idxC = 1: cMvCell
                		energyMv = energyMv + sqrt(sum(mvCell{idxR, idxC} .^ 2));
                	end
                end
                energyMv = energyMv / rMvCell / cMvCell;
                sumEnergyMvAllLibPic(idxLibPic) = sumEnergyMvAllLibPic(idxLibPic) + energyMv;
			end

            averTextureStd = averTextureStd / numIntraPic;

			sumFreq = sum(listLibFreq - 1); % remove the library picture itself.
	        weight = (listLibFreq - 1) / sumFreq;
            weight = weight';
	        expectCost = sum(weight .* listClusterCost);
	        expectDist = sum(weight .* sumDistAllLibPic);
            expectSim = sum(weight .* sumSimAllLibPic);
            expectDistPoc = sum(weight .* sumDistPocAllLibPic);
            expectSimPoc = sum(weight .* sumSimPocAllLibPic);
            expectDistMC = sum(weight .* sumDistMCAllLibPic);
            expectEnergyMv = sum(weight .* sumEnergyMvAllLibPic);

            contentInfoAllSeq(idxSeq, 1, 1) = averTextureStd;            
			contentInfoAllSeq(idxSeq, 1, 2) = expectCost;
			contentInfoAllSeq(idxSeq, 2, 2) = expectDist;
            contentInfoAllSeq(idxSeq, 3, 2) = expectSim;
            contentInfoAllSeq(idxSeq, 4, 2) = expectDistPoc;
            contentInfoAllSeq(idxSeq, 5, 2) = expectSimPoc;
            contentInfoAllSeq(idxSeq, 6, 2) = expectDistMC;
            contentInfoAllSeq(idxSeq, 7, 2) = expectEnergyMv;
		end
    end 

    %% compute OPT-DQP for samples.
	if subFlagBestQP
		numBestDqp = 3;
		[gainSort, idxSort] = sort(gainYUVAllSeq, 2, 'ascend');
		bestDqpSeqAllSeq = zeros(numAllSeq, numBestDqp);
        weightBestDqpSeqAllSeq = zeros(numAllSeq, numBestDqp);
        rangeWeight = 0.005;
		for idxSeq = 1: numAllSeq
			idxBestDqp = idxSort(idxSeq, 1: numBestDqp);
			bestDqpSeqAllSeq(idxSeq, :) = listDqp(idxBestDqp);
			bestPerformance = gainYUVAllSeq(idxSeq, idxBestDqp);
            weightBestDqpSeqAllSeq(idxSeq, :) = max(0, 1 - (bestPerformance - bestPerformance(1)) / rangeWeight);
        end
        
        bestDqpAllSeq = bestDqpSeqAllSeq;
        weightBestDqpAllSeq = weightBestDqpSeqAllSeq;
    end
    
    if subFlagBestQPKeyFrame
		numBestDqp = 3;
		[gainSort, idxSort] = sort(gainKeyYUVAllSeq, 2, 'ascend');
		bestDqpKeyAllSeq = zeros(numAllSeq, numBestDqp);
        weightBestDqpKeyAllSeq = zeros(numAllSeq, numBestDqp);
        rangeWeight = 0.01;
		for idxSeq = 1: numAllSeq
			idxBestDqp = idxSort(idxSeq, 1: numBestDqp);
			bestDqpKeyAllSeq(idxSeq, :) = listDqp(idxBestDqp);
			bestPerformance = gainKeyYUVAllSeq(idxSeq, idxBestDqp);
            weightBestDqpKeyAllSeq(idxSeq, :) = max(0, 1 - (bestPerformance - bestPerformance(1)) / rangeWeight);
        end
        bestDqpAllSeq = bestDqpKeyAllSeq;
        weightBestDqpAllSeq = weightBestDqpKeyAllSeq;
    end
    

    %% compute affecting factor for samples.
	if subFlagShowData
		subsubFlag1 = 1;
		subsubFlag2 = 1;
        subsubFlag3 = 1;
		subsubFlag4 = 1;
		subsubFlag5 = 1;
        subsubFlag6 = 1;

		sampleSet = zeros(numAllSeq * numBestDqp, 3);
		weightSet = zeros(numAllSeq * numBestDqp, 1);
		
		% 1----average std and expectation cluster cost.
		if subsubFlag1
			for idxBestDqp = 1: numBestDqp
				startPos = (idxBestDqp - 1) * numAllSeq + 1;
				endPos = startPos + numAllSeq - 1;
				sampleSet(startPos: endPos, :) = [contentInfoAllSeq(:, 1, 1) contentInfoAllSeq(:, 1, 2) bestDqpAllSeq(:, idxBestDqp)];
				weightSet(startPos: endPos) = weightBestDqpAllSeq(:, idxBestDqp);
            end
            x1 = sampleSet(:, 1);
            y1 = sampleSet(:, 2);
            z1 = sampleSet(:, 3);
            x1_1 = log(sampleSet(:, 1) + 1);
            y1_1 = log(sampleSet(:, 2) + 1);
            z1_1 = sampleSet(:, 3);
        end
        % 2----average std and expectation frame difference.
		if subsubFlag2
			for idxBestDqp = 1: numBestDqp
				startPos = (idxBestDqp - 1) * numAllSeq + 1;
				endPos = startPos + numAllSeq - 1;
				sampleSet(startPos: endPos, :) = [contentInfoAllSeq(:, 1, 1) contentInfoAllSeq(:, 2, 2) bestDqpAllSeq(:, idxBestDqp)];
            	weightSet(startPos: endPos) = weightBestDqpAllSeq(:, idxBestDqp);
            end
            x2 = sampleSet(:, 1);
            y2 = sampleSet(:, 2);
            z2 = sampleSet(:, 3);
            x2_1 = log(sampleSet(:, 1) + 1);
            y2_1 = log(sampleSet(:, 2) + 1);
            z2_1 = sampleSet(:, 3);
        end
		% 3----average std and expectation similarity based on frame difference.
		if subsubFlag3
			for idxBestDqp = 1: numBestDqp
				startPos = (idxBestDqp - 1) * numAllSeq + 1;
				endPos = startPos + numAllSeq - 1;
				sampleSet(startPos: endPos, :) = [contentInfoAllSeq(:, 1, 1) contentInfoAllSeq(:, 3, 2) bestDqpAllSeq(:, idxBestDqp)];
            	weightSet(startPos: endPos) = weightBestDqpAllSeq(:, idxBestDqp);
            end
            x3 = sampleSet(:, 1);
            y3 = sampleSet(:, 2);
            z3 = sampleSet(:, 3);
            x3_1 = log(sampleSet(:, 1) + 1);
            y3_1 = log(sampleSet(:, 2) + 1);
            z3_1 = sampleSet(:, 3);
        end
        % 4----average std and expect frame poc difference.
		if subsubFlag4
			for idxBestDqp = 1: numBestDqp
				startPos = (idxBestDqp - 1) * numAllSeq + 1;
				endPos = startPos + numAllSeq - 1;
				sampleSet(startPos: endPos, :) = [contentInfoAllSeq(:, 1, 1) contentInfoAllSeq(:, 4, 2) bestDqpAllSeq(:, idxBestDqp)];
            	weightSet(startPos: endPos) = weightBestDqpAllSeq(:, idxBestDqp);
            end
            x4 = sampleSet(:, 1);
            y4 = sampleSet(:, 2);
            z4 = sampleSet(:, 3);
            x4_1 = log(sampleSet(:, 1) + 1);
            y4_1 = log(sampleSet(:, 2) + 1);
            z4_1 = sampleSet(:, 3);
        end
        % 5----average std and expect similarity based on frame poc difference.
		if subsubFlag5
			for idxBestDqp = 1: numBestDqp
				startPos = (idxBestDqp - 1) * numAllSeq + 1;
				endPos = startPos + numAllSeq - 1;
				sampleSet(startPos: endPos, :) = [contentInfoAllSeq(:, 1, 1) contentInfoAllSeq(:, 5, 2) bestDqpAllSeq(:, idxBestDqp)];
            	weightSet(startPos: endPos) = weightBestDqpAllSeq(:, idxBestDqp);
            end
            x5 = sampleSet(:, 1);
            y5 = sampleSet(:, 2);
            z5 = sampleSet(:, 3);
            x5_1 = log(sampleSet(:, 1) + 1);
            y5_1 = log(sampleSet(:, 2) + 1);
            z5_1 = sampleSet(:, 3);
        end
        % 6----average std and expect difference based on motion compensation.
		if subsubFlag6
			for idxBestDqp = 1: numBestDqp
				startPos = (idxBestDqp - 1) * numAllSeq + 1;
				endPos = startPos + numAllSeq - 1;
				sampleSet(startPos: endPos, :) = [contentInfoAllSeq(:, 1, 1) contentInfoAllSeq(:, 6, 2) bestDqpAllSeq(:, idxBestDqp)];
            	weightSet(startPos: endPos) = weightBestDqpAllSeq(:, idxBestDqp);
            end
            x6 = sampleSet(:, 1);
            y6 = sampleSet(:, 2);
            z6 = sampleSet(:, 3);
            x6_1 = log(sampleSet(:, 1) + 1);
            y6_1 = log(sampleSet(:, 2) + 1);
            z6_1 = sampleSet(:, 3);
        end
	end
end