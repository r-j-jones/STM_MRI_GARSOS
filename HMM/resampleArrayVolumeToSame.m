function resampledVolume = resampleArrayVolumeToSame(srcVolume, refVolume, varargin)

nRefDimensions = max(4, ndims(refVolume));
refArraySize = cell(1, nRefDimensions);
[refArraySize{:}] = size(refVolume.A);
refArraySize = cell2mat(refArraySize);
    
nSrcDimensions = max(4, ndims(srcVolume));
srcArraySize = cell(1, nSrcDimensions);
[srcArraySize{:}] = size(srcVolume.A);
srcArraySize = cell2mat(srcArraySize);

% assert(all(refArraySize(4:end) == srcArraySize(4:end)), 'Non-spatial dimensions must have the same size.');

ii = arrayfun(@(x)0:x-1, refArraySize(1:3),'UniformOutput',false);
II = cell(3, 1);
[II{:}] = ndgrid(ii{:});
II = cellfun(@(x)reshape(x,1,[]),II,'UniformOutput',false);
II = cell2mat(II);

JJ = bsxfun(@rdivide,srcVolume.R\bsxfun(@plus,refVolume.R*bsxfun(@times,refVolume.v,II),refVolume.r0 - srcVolume.r0),srcVolume.v) + 1 + 1;

JJ = mat2cell(JJ,[1 1 1],prod(refArraySize(1:3)));

JJ = cellfun(@(x)reshape(x,refArraySize(1:3)), JJ, 'UniformOutput', false);

resArraySize = [refArraySize(1:3) srcArraySize(4:end)];

B = zeros(resArraySize);

for iBlock = 1:prod(resArraySize(4:end))
    B(:,:,:,iBlock) = interpn(centeredResize(srcVolume.A(:,:,:,iBlock),srcArraySize(1:3) + 2), JJ{:}, varargin{:});
end

resampledVolume = ArrayVolume(B, refVolume.R, refVolume.v, refVolume.r0, []);






























