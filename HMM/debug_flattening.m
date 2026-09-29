rawKxy = rawLocations4(:, :, :, 1:2);
transformedKxy = locations4(:, :, :, 1:2);

rawPartitionDifference = max(abs( ...
    rawKxy - repmat(rawKxy(:, 1, :, :), [1 nPartition 1 1])), [], 'all');

transformedPartitionDifference = max(abs( ...
    transformedKxy - repmat(transformedKxy(:, 1, :, :), ...
    [1 nPartition 1 1])), [], 'all');

fprintf('Raw partition Kxy difference: %.17g\n', ...
    rawPartitionDifference);
fprintf('Transformed partition Kxy difference: %.17g\n', ...
    transformedPartitionDifference);

fprintf('CA =\n');
disp(coordinateTransform);

fprintf('CA partition-to-in-plane terms: [% .17g % .17g]\n', ...
    coordinateTransform(3, 1), coordinateTransform(3, 2));

fprintf('nReadout=%d, nPartition=%d, nSpokes=%d\n', ...
    nReadout, nPartition, nSpokes);






disp(squeeze(rawLocations4(1, 1:min(nPartition,5), 1, :)));
disp(squeeze(rawLocations4(1, 1:min(nPartition,5), 2, :)));
disp(squeeze(locations4(1, 1:min(nPartition,5), 1, :)));
disp(squeeze(locations4(1, 1:min(nPartition,5), 2, :)));


% disp(squeeze(rawLocations4(end, 1:min(nPartition,5), 1, :)));
% disp(squeeze(rawLocations4(end, 1:min(nPartition,5), 2, :)));
% disp(squeeze(locations4(end, 1:min(nPartition,5), 1, :)));
% disp(squeeze(locations4(end, 1:min(nPartition,5), 2, :)));


difference = abs(locations4(:, :, :, 1:2) - expectedKxy);
[maxDifference, linearIndex] = max(difference(:));

fprintf('Maximum Kxy mismatch: %.17g\n', maxDifference);
fprintf('Mismatch linear index: %d\n', linearIndex);


