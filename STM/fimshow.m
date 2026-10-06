function fig = fimshow( imdata, colormapToUse, addColorbar, titleString )

if nargin<2
    colormapToUse = [];
end
if nargin<3
    addColorbar = false;
end
if nargin<4
    titleString = [];
end

 % Create a figure
fig = figure();

% Display image and fit it to the window
imshow((imdata), [], 'InitialMagnification', 'fit');

% Remove outer padding and make axes fill the whole figure
ax = gca;
ax.Units = 'normalized';
ax.Position = [0 0 1 1]; 

if ~isempty(colormapToUse)
    colormap(colormapToUse);
end

if addColorbar 
    ax.Position = [0 0.05 0.9 0.9];
    cbb = colorbar;
    cbb.FontSize = 13;
end

if ~isempty(titleString)
    if ax.Position(2) > 0
        ax.Position(2) = 0.01; %ax.Position(2)-0.05;
    end
    title(titleString,'FontSize',18);
end

% Optional: Hide axis ticks/borders completely
axis off;

drawnow;

end
