function [err] = diffBlock(refBlock, trgBlock)

% average difference.
err = sum((refBlock(:) - trgBlock(:)).^2);

end