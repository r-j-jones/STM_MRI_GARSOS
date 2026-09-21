function X = fftlattice(A, b, n)

xi = arrayfun(@(x){(1:x) - (floor(x/2) + 1)}, n);

X = cell(1, numel(b));
[X{:}] = ndgrid(xi{:});

X = cellfun(@(x){reshape(x,[1 prod(n)])}, X);

X = cat(1, X{:});

X = bsxfun(@plus, A*X, b);
