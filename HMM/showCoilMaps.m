function f = showCoilMaps(data,dataDescription,nDimensions)

inputArraySize = ones(1,max(nDimensions + 1,ndims(data)));
inputArraySize(1:ndims(data)) = size(data);
inputArraySize = [inputArraySize(1,1:nDimensions) prod(inputArraySize)/prod(inputArraySize(1,1:nDimensions))];

nCoils = inputArraySize(end);

arraySize = ones(1,max(3,nDimensions));
arraySize(1:nDimensions) = inputArraySize(1:end-1);

arraySize = [arraySize(1,1:2) prod(arraySize)/prod(arraySize(1:2))];

data = reshape(data,[arraySize nCoils]);

iSlice = floor(arraySize(3)/2) + 1;

f = figure('NumberTitle', 'off','name',dataDescription); %,'MenuBar','none','Toolbar','none');
set(f,'Color',[0.5 0.5 0.5]);

imageHandles = cell(1, nCoils);
axesHandles = cell(1, nCoils);

for iCoil = 1:nCoils
    %indices = arrayfun(@(x,y)x:y, [ones(1,numel(arraySize)), i], [arraySize, i],'UniformOutput',false);
    axesHandles{iCoil} = subplot(ceil(nCoils/ceil(sqrt(nCoils))),ceil(sqrt(nCoils)),iCoil);
    imageHandles{iCoil} = imagesc(complexToColor(data(:,:,iSlice,iCoil),'log',true));
    set(axesHandles{iCoil},'XColor',[0.5 0.5 0.5]*1.5);
    set(axesHandles{iCoil},'YColor',[0.5 0.5 0.5]*1.5);
    title(sprintf('Slice %.0f',iSlice),'Color',[0.5 0.5 0.5]*1.5);
    axis equal; axis tight; axis off;
    %axis tight;
end

%maximizeFigure;

figureGuiData = guidata(f);
figureGuiData.data = data;
figureGuiData.nCoils = nCoils;
figureGuiData.arraySize = arraySize;
figureGuiData.iSlice = iSlice;
figureGuiData.imageHandles = imageHandles;
figureGuiData.axesHandles = axesHandles;
guidata(f,figureGuiData);

set(f,'WindowScrollWheelFcn',@showCoilImagesScrollFunction);

end

function showCoilImagesScrollFunction(hObject, eventdata, handles, varargin) %#ok

figureGuiData = guidata(hObject);

figureGuiData.iSlice = max(1,min(figureGuiData.iSlice + eventdata.VerticalScrollCount,figureGuiData.arraySize(3)));
guidata(hObject,figureGuiData);

for iCoil = 1:figureGuiData.nCoils
    set(figureGuiData.imageHandles{iCoil},'CData',complexToColor(figureGuiData.data(:,:,figureGuiData.iSlice,iCoil),'log',true));
    set(get(figureGuiData.axesHandles{iCoil},'Title'),'String',sprintf('Slice %.0f',figureGuiData.iSlice));
end

end



















