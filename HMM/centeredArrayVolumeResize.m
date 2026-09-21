function b = centeredArrayVolumeResize(a, newArraySize)

assert(all(size(newArraySize) == [1, 3]));

nDimensions = max(ndims(a.A), 3);

oldArraySize = cell(1, nDimensions);
[oldArraySize{:}] = size(a.A);

oldArraySize = cell2mat(oldArraySize);

newArraySize = [newArraySize, oldArraySize(1,4:end)];

B = centeredResize(a.A, newArraySize);

cOld = floor(oldArraySize(1,1:3)/2) + 1;
cNew = floor(newArraySize(1,1:3)/2) + 1;

r0 = a.r0 + a.R*diag(a.v)*(cOld - cNew)';

b = ArrayVolume(B, a.R, a.v, r0, []);

