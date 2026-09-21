function arraySize = getArraySize(A, n)

arraySize = cell(1, n);

[arraySize{:}] = size(A);

arraySize = cell2mat(arraySize);