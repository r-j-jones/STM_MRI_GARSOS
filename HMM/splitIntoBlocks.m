function blocks = splitIntoBlocks(a, m)

n = numel(a);

iiStart = 1:m:n;
iiEnd = iiStart + m - 1;
iiEnd = min(iiEnd, n);

blocks = arrayfun(@(x, y)a(1, x:y), iiStart, iiEnd, 'UniformOutput', false);