

%% TEST 1

optest = result.stmOperator.radialOperator;

xTest = complex(randn(224,224,133,'single'),...
    randn(224,224,133,'single'));

t = 1;

coilImages = reshape( ...
    xTest(:, :, t) .* optest.senseMaps, ...
    prod(optest.imageSize), optest.nCoil);


yBatch = nufft(coilImages, optest.plans(t).st);
    disp(size(yBatch));


coilImages = xTest(:, :, t) .* optest.senseMaps;
size(coilImages)
% Expected: 224 x 224 x 20

yBatch = nufft(coilImages, optest.plans(t).st);
size(yBatch)


yMatrix = reshape(yBatch, ...
    optest.nReadout * optest.nSpokes, optest.nCoil);

xBatch = nufft_adj(yMatrix, optest.plans(t).st);
size(xBatch)


%% TEST 2

coilImages = complex( ...
    randn(224,224,20,'single'), ...
    randn(224,224,20,'single'));

y = nufft(coilImages, result.radialOperator.plans(1).st);

fprintf('Input class: %s\n', class(coilImages));
fprintf('sn class: %s\n', class(result.radialOperator.plans(1).st.sn));
fprintf('Output class: %s\n', class(y));



stSingle = result.radialOperator.plans(1).st;
stSingle.sn = single(stSingle.sn);

ySingleSn = nufft(coilImages, stSingle);

difference = norm(double(y(:)) - double(ySingleSn(:))) / ...
    max(norm(double(y(:))), eps);

fprintf('Relative difference: %.3e\n', difference);
fprintf('Output with single sn: %s\n', class(ySingleSn));



%% BAD TEST


% yTest = optest.forward(xTest);
% 
% whos xTest yTest
% class(optest.plans(1).omega)
% 
% t = 1;
% 
% coilImages = reshape( ...
%     xTest(:, :, t) .* senseMaps, ...
%     prod(imageSize), nCoil);
% 
% try
%     yBatch = nufft(coilImages, plans(t).st);
%     disp(size(yBatch));
% catch ME
%     disp(getReport(ME, 'extended'));
% end