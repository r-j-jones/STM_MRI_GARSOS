function check_dcf( densityCompensation, nPartitions, nSpokes, nSamples)

dcf3 = reshape(densityCompensation,[nSamples nPartitions nSpokes]);

reference = dcf3(:,:,1); 

%[[[[ Check 1
maxDiff = max(abs(dcf3 - reference), [], 'all'); 
fprintf('Maximum absolute difference: %.3e\n', maxDiff);

%[[[[ Check 2
maxRelDiff = max(abs(dcf3 - reference) ./ max(abs(reference), eps), [], 'all');
fprintf('Maximum relative difference: %.3e\n', maxRelDiff);

%[[[[ Check 3
dcf3 = reshape(densityCompensation, [384, 58, 2000]);

% Variation over dimension 1
varDim1 = max(abs(dcf3 - dcf3(1,:,:)), [], 'all');

% Variation over dimension 2
varDim2 = max(abs(dcf3 - dcf3(:,1,:)), [], 'all');

% Variation over dimension 3
varDim3 = max(abs(dcf3 - dcf3(:,:,1)), [], 'all');

fprintf('Variation relative to index 1:\n');
fprintf('  dim 1: %.3e\n', varDim1);
fprintf('  dim 2: %.3e\n', varDim2);
fprintf('  dim 3: %.3e\n', varDim3);

%[[[[ Check 4
reference = dcf3(:,:,1);

differences = abs(dcf3 - reference);

maxDifferencePerSpoke = squeeze(max(max(differences, [], 1), [], 2));

fprintf('Maximum difference for each spoke:\n');
disp(maxDifferencePerSpoke.');


%[[[[ Check 5
dcf3 = reshape(densityCompensation, [384, 58, 2000]);

reference = dcf3(:,:,1);

maxDiffPerSpoke = zeros(2000,1);

for i = 1:2000
    maxDiffPerSpoke(i) = max(abs(dcf3(:,:,i) - reference), [], 'all');
end

fprintf('Maximum difference across all spokes: %.6e\n', max(maxDiffPerSpoke));
fprintf('Number of exactly identical spokes: %d / %d\n', ...
    sum(maxDiffPerSpoke == 0), 2000);

