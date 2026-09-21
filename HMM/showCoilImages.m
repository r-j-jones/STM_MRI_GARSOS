function f = showCoilImages(data,dataDescription,nDimensions)

inputArraySize = ones(1, max(nDimensions + 1, ndims(data)));
inputArraySize(1:ndims(data)) = size(data);
inputArraySize = [inputArraySize(1,1:nDimensions) prod(inputArraySize)/prod(inputArraySize(1,1:nDimensions))];

nCoils = inputArraySize(end);

arraySize = ones(1,max(3,nDimensions));
arraySize(1:nDimensions) = inputArraySize(1:end-1);

arraySize = [arraySize(1,1:2) prod(arraySize)/prod(arraySize(1:2))];

data = reshape(data,[arraySize nCoils]);

iSlice = floor(arraySize(3)/2) + 1;

f = figure('NumberTitle', 'off', 'name', dataDescription, 'Color', 'k');

imageHandles = cell(1, nCoils);
axesHandles = cell(1, nCoils);

maxValue = max(abs(data(:)));

for iCoil = 1:nCoils
    %indices = arrayfun(@(x,y)x:y, [ones(1,numel(arraySize)), i], [arraySize, i],'UniformOutput',false);
    axesHandles{iCoil} = subplot(ceil(nCoils/ceil(sqrt(nCoils))),ceil(sqrt(nCoils)),iCoil);
    imageHandles{iCoil} = image(complexToColor(data(:,:,iSlice,iCoil), 'max', maxValue));
    set(axesHandles{iCoil},'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5]);
    title(sprintf('Channel %.0f, Slice %.0f',iCoil,iSlice),'Color',[0.5 0.5 0.5]);
    axis equal; axis tight;
%     axis normal; axis tight;
end

%maximizeFigure;

figureGuiData = guidata(f);
figureGuiData.data = data;
figureGuiData.nCoils = nCoils;
figureGuiData.arraySize = arraySize;
figureGuiData.iSlice = iSlice;
figureGuiData.imageHandles = imageHandles;
figureGuiData.axesHandles = axesHandles;
figureGuiData.ctrlIsDown = false;
figureGuiData.intensityScale = 0;
figureGuiData.maxValue = maxValue;
guidata(f,figureGuiData);

set(f,'WindowScrollWheelFcn',@showCoilImagesScrollFunction);
set(f,'WindowKeyPressFcn',@showCoilImagesKeyPressFunction);
set(f,'WindowKeyReleaseFcn',@showCoilImagesKeyReleaseFunction);

for iCoil = 1:nCoils
    set(imageHandles{iCoil},'ButtonDownFcn',@showCoilImagesButtonDownFunction);
end

end

function showCoilImagesScrollFunction(hObject, eventdata, handles, varargin) 

figureGuiData = guidata(hObject);

if figureGuiData.ctrlIsDown
%     fprintf('Change slice orientation.\n');
    figureGuiData.intensityScale = figureGuiData.intensityScale - eventdata.VerticalScrollCount;
else
%     fprintf('Change slice position.\n');
    figureGuiData.iSlice = max(1,min(figureGuiData.iSlice - eventdata.VerticalScrollCount,figureGuiData.arraySize(3)));
end

guidata(hObject,figureGuiData);

for iCoil = 1:figureGuiData.nCoils
    A = complexToColor(figureGuiData.data(:,:,figureGuiData.iSlice,iCoil),'max',figureGuiData.maxValue*1.05^-figureGuiData.intensityScale);
    
%     A = A*(1.05^figureGuiData.intensityScale);
%     
%     A(A > 1) = 1;
%     
%     A(:,:,1) = A(:,:,2);
%     A(:,:,3) = A(:,:,2);

%     imwrite(A,fullfile('c:\adamroot\tmp\test102',sprintf('Instance%.0fCoil%.0fSlice%.0f.png',1000*gcf,iCoil,figureGuiData.iSlice)));
    %fprintf('Olle');
    set(figureGuiData.imageHandles{iCoil},'CData',A);
    set(get(figureGuiData.axesHandles{iCoil},'Title'),'String',sprintf('Channel %.0f, Slice %.0f',iCoil,figureGuiData.iSlice));
end

end

function showCoilImagesButtonDownFunction(hObject, eventdata, handles, varargin)

figureGuiData = guidata(hObject);

for iCoil = 1:numel(figureGuiData.imageHandles)
    
    if hObject == figureGuiData.imageHandles{iCoil}
        
        p = get(figureGuiData.axesHandles{iCoil},'CurrentPoint');
        
        x = round(p(1,1));
        y = round(p(1,2));
        
        value = figureGuiData.data(y,x,figureGuiData.iSlice,iCoil);
        
%         fprintf('Clicked! (%.0f,%.0f)\n',x,y);
        
        set(get(figureGuiData.axesHandles{iCoil},'Title'),'String',sprintf('Channel %.0f, Slice %.0f (%s)',iCoil,figureGuiData.iSlice,num2str(value)));
        
    end
    
end

end

function showCoilImagesKeyPressFunction(hObject, eventdata, handles, varargin)
    figureGuiData = guidata(hObject);
    figureGuiData.ctrlIsDown = true;
    guidata(hObject,figureGuiData);
end

function showCoilImagesKeyReleaseFunction(hObject, eventdata, handles, varargin)
    figureGuiData = guidata(hObject);
    figureGuiData.ctrlIsDown = false;
    guidata(hObject,figureGuiData);
end



















