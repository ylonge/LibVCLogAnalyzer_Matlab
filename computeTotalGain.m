function [ gainAllSeq ] = computeTotalGain(bitPsnrEachLibALLSeqLibvc, bitPsnrALLSeqAnchor, numFrameEachLibAllSeq, listDqp)
%Author: ylonge.
%Function: compute fit parameter and show curve with given equation. Only one independent factor is considered.
%-Input:
%   --bitPsnrEachLibALLSeqLibvc: storing all bit psnr of frames referencing each library picture.
%	--bitPsnrALLSeqAnchor: bit psnr of anchor.
%	--numFrameEachLibAllSeq: frame number of key pictures that referencing one library picture.
%	--listDqp: the given delta QP.
%-OutPut:
%	--gainAllSeq: the total gain.

%% computation.
% prepare data.
numAllSeq = size(bitPsnrEachLibALLSeqLibvc, 1);
countTemp = 0;

% collect data.
gainAllSeq = zeros(numAllSeq, 4);

% for each sequence.
for idxSeq = 1: numAllSeq
	% prepare data.
	bitPsnrAnchor = bitPsnrALLSeqAnchor{idxSeq};
	numFrameEachLib = numFrameEachLibAllSeq{idxSeq};
	numLibPic = length(numFrameEachLib);

	% sum bit psnr of libvc.
	bitPsnrLibvc = zeros(4, 5);
    for idxLibPic = 1: numLibPic
    	countTemp = countTemp + 1;
    	idxDqp = 3 - listDqp(countTemp); % for dqp= 2~-11;
    	interval = ((idxLibPic - 1)*4 + 1): (idxLibPic * 4);
    	bitPsnrEachLib = bitPsnrEachLibALLSeqLibvc{idxSeq, idxDqp}(interval, :);
    	bitPsnrLibvc = bitPsnrLibvc + bitPsnrEachLib * numFrameEachLib(idxLibPic);
    end
    bitPsnrLibvc = bitPsnrLibvc / sum(numFrameEachLib);

    % BD-rate.
    gainAllSeq(idxSeq, :) = bdRateComparation(bitPsnrAnchor, bitPsnrLibvc);
end

end