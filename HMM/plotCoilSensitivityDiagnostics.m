function figs = plotCoilSensitivityDiagnostics(debugStructure, sliceNumbers)
%PLOTCOILSENSITIVITYDIAGNOSTICS Display orthogonal diagnostic slices.
%
% Creates separate figures for the surface-coil combined image and the
% body-coil reference image.
%
% ARRAY DIMENSION ORDER
%   The volume is stored as volume(Y,Z,X):
%       dimension 1 = Y
%       dimension 2 = Z
%       dimension 3 = X
%
% Each figure uses a horizontal 1-by-3 subplot arrangement:
%       X-normal slice | Y-normal slice | Z-normal slice
%
% INPUTS
%   debugStructure
%       Structure containing T2 and bodyCoilImage.
%
%   sliceNumbers (optional)
%       Three positive integers [xSlice,ySlice,zSlice]. Center slices are
%       used by default.
%
% OUTPUT
%   figs
%       Structure containing the generated figure handles.

arguments
    debugStructure (1,1) struct
    sliceNumbers (1,3) double {mustBeInteger, mustBePositive} = []
end

fieldNames = {'T2', 'bodyCoilImage'};

descriptions = { ...
    'Surface-coil combined image', ...
    'Body-coil reference image'};

figs = struct();

for volumeIndex = 1:numel(fieldNames)

    fieldName = fieldNames{volumeIndex};
    figureName = descriptions{volumeIndex};

    if ~isfield(debugStructure, fieldName)
        error('debugStructure is missing the required field "%s".', ...
            fieldName);
    end

    % Extract the numerical data from the ArrayVolume object.
    volume = squeeze(debugStructure.(fieldName).A);

    if ndims(volume) ~= 3
        error('The diagnostic volume "%s" must be three-dimensional.', ...
            fieldName);
    end

    volumeSize = size(volume);

    % Array dimension mapping:
    %   dimension 1 = Y
    %   dimension 2 = Z
    %   dimension 3 = X
    xDimension = 3;
    yDimension = 1;
    zDimension = 2;

    if isempty(sliceNumbers)
        xSlice = floor(volumeSize(xDimension) / 2) + 1;
        ySlice = floor(volumeSize(yDimension) / 2) + 1;
        zSlice = floor(volumeSize(zDimension) / 2) + 1;
    else
        xSlice = sliceNumbers(1);
        ySlice = sliceNumbers(2);
        zSlice = sliceNumbers(3);
    end

    if xSlice > volumeSize(xDimension) || ...
       ySlice > volumeSize(yDimension) || ...
       zSlice > volumeSize(zDimension)

        error([ ...
            'Requested [X Y Z] slices [%d %d %d] exceed the available ' ...
            '[X Y Z] dimensions [%d %d %d] for "%s".'], ...
            xSlice, ySlice, zSlice, ...
            volumeSize(xDimension), ...
            volumeSize(yDimension), ...
            volumeSize(zDimension), ...
            fieldName);
    end

    % Display magnitude for both real and complex volumes.
    displayVolume = abs(volume);

    % Extract the three displayed slices before plotting so that common
    % colormap limits can be calculated across all three panels.
    xSliceData = squeeze(displayVolume(:, :, xSlice))';
    ySliceData = squeeze(displayVolume(ySlice, :, :));
    zSliceData = squeeze(displayVolume(:, zSlice, :));

    % Use the minimum and maximum values across all displayed slices.
    commonColorLimits = [ ...
        min([xSliceData(:); ySliceData(:); zSliceData(:)]), ...
        max([xSliceData(:); ySliceData(:); zSliceData(:)])];

    % Avoid invalid identical color limits for constant-valued data.
    if commonColorLimits(1) == commonColorLimits(2)
        commonColorLimits(2) = commonColorLimits(1) + eps( ...
            commonColorLimits(1));
    end

    fig = figure( ...
        'Name', figureName, ...
        'Color', 'w', ...
        'Position', [385 366 1388 481]);

    layout = tiledlayout(fig, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');


    %% X-normal slice: fix dimension 3 and display Y versus Z

    axX = nexttile(layout);

    imagesc(axX, xSliceData, commonColorLimits);

    axis(axX, 'image');
    axis(axX, 'xy');
    axX.FontSize = 12;

    colorbar(axX, 'FontSize', 12);

    xlabel(axX, 'Y index', 'FontSize', 13);
    ylabel(axX, 'Z index', 'FontSize', 13);

    title(axX, sprintf('X slice %d', xSlice), 'FontSize', 14);


    %% Y-normal slice: fix dimension 1 and display X versus Z

    axY = nexttile(layout);

    imagesc(axY, ySliceData, commonColorLimits);

    axis(axY, 'image');
    axis(axY, 'xy');
    axY.FontSize = 12;

    colorbar(axY, 'FontSize', 12);

    xlabel(axY, 'X index', 'FontSize', 13);
    ylabel(axY, 'Z index', 'FontSize', 13);

    title(axY, sprintf('Y slice %d', ySlice), 'FontSize', 14);


    %% Z-normal slice: fix dimension 2 and display X versus Y

    axZ = nexttile(layout);

    imagesc(axZ, zSliceData, commonColorLimits);

    axis(axZ, 'image');
    axis(axZ, 'xy');
    axZ.FontSize = 12;

    colorbar(axZ, 'FontSize', 12);

    xlabel(axZ, 'X index', 'FontSize', 13);
    ylabel(axZ, 'Y index', 'FontSize', 13);

    title(axZ, sprintf('Z slice %d', zSlice), 'FontSize', 14);


    %% Overall title and colormap

    title(layout, figureName, 'FontSize', 18);

    if exist('viridis', 'file') == 2
        colormap(fig, viridis(256));
    else
        colormap(fig, parula(256));
    end

    figs.(fieldName) = fig;
end

end