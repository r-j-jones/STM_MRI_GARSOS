function variableValues = getIndexVariablesFromDicomHeaders(dicomHeaders, variableNames)

t = cellfun(@(x)getIndexVariablesFromDicomHeader(x, variableNames), dicomHeaders,'UniformOutput',false);

variableValues = cat(1, t{:});

function variableValues = getIndexVariablesFromDicomHeader(dicomHeader, variableNames)

nVariables = numel(variableNames);

variableValues = cell(1, nVariables);

for iVariable = 1:nVariables
    switch variableNames{iVariable}
        case 'B_VALUE'            
            %bValue = sscanf(char(dicomHeader.Private_0019_100c'), '%d');
            
            bValue = sprintf('%012d', dicomHeader.SiemensBValue);
            
            %variableValues{iVariable} = sscanf(char(dicomHeader.Private_0019_100c'),'%d');
            
            diffusionDirectionality = dicomHeader.SiemensDiffusionDirectionality;
            
            switch diffusionDirectionality
                case {'NONE','ISOTROPIC'}
                    diffusionGradientDirectionIndexStub = diffusionDirectionality;
                case 'DIRECTIONAL'
                    diffusionGradientDirection = typecast(dicomHeader.Private_0019_100e','double');
                    diffusionGradientDirectionIndexStub = sprintf('%s%s',diffusionDirectionality,sprintf('_%.5f', diffusionGradientDirection));
                otherwise
                    error('Unsupported diffusion directionality.');
            end
            
            variableValues{iVariable} = {sprintf('%s_%s',bValue, diffusionGradientDirectionIndexStub)};
            
        case 'GE_B_VALUE'
            
            geBValueString = char(dicomHeader.Private_0043_1039');
            
            geBValue = sscanf(geBValueString,'%f');
            
            bValue = sprintf('%012d', geBValue);
            
            variableValues{iVariable} = {sprintf('%s_%s',bValue, 'GE_B_VALUE')};
            
        case 'GE_B_VALUE_INSTANCE_NUMBER'
            
            %geBValueInstanceNumber = floor((dicomHeader.InstanceNumber - 1)/46);
            
            %bValueInstanceNumber = sprintf('%012d', dicomHeader.InstanceNumber);
            
            %variableValues{iVariable} = {sprintf('%s_%s', bValueInstanceNumber, 'GE_B_VALUE_INSTANCE_NUMBER')};
            
            variableValues{iVariable} = dicomHeader.InstanceNumber;
            
        case 'PHILIPS_B_VALUE'
            
        
        case 'AcquisitionNumber'
            variableValues{iVariable} = dicomHeader.AcquisitionNumber;
        case 'InstanceNumber'
            variableValues{iVariable} = dicomHeader.InstanceNumber;
        case 'DiffusionBValue'
            variableValues{iVariable} = dicomHeader.DiffusionBValue;
        case 'TTC'
            t = sscanf(dicomHeader.ImageComments, 'TTC %e sec');            
            variableValues{iVariable} = {sprintf('TTC %025.10f sec', t)};
        case 'RADIAL_VIBE_SERIES_SUFFIX'
            [~, radialVibeSeriesSuffixNumber] = regexp(dicomHeader.SeriesDescription, '.+_T([0-9]+)\s*$','match','tokens','once');
            radialVibeSeriesSuffixNumber = sscanf(radialVibeSeriesSuffixNumber{1}, '%e');
            variableValues{iVariable} = {sprintf('RADIAL_VIBE_SERIES_SUFFIX %08.0f', radialVibeSeriesSuffixNumber)};
        case 'RADIAL_VIBE_ACQUISITION_TIME_VE11'
            variableValues{iVariable} = {sprintf('RADIAL_VIBE_SERIES_SUFFIX_VE11 %025.10f s', convertAcquisitionTimeToSeconds(dicomHeader.AcquisitionTime))};
        case 'RADIAL_VIBE_SERIES_SUFFIX_VE11'
            %[~, radialVibeSeriesSuffixNumber] = regexp(dicomHeader.SeriesDescription, '.+_TT=([0-9.]+)s.*$','match','tokens','once');
            %[~, radialVibeSeriesSuffixNumber] = regexp(dicomHeader.SeriesDescription, '.+_TT=([0-15.]+)s.*$','match','tokens','once');
            %radialVibeSeriesSuffixNumber = sscanf(radialVibeSeriesSuffixNumber{1}, '%e');
            radialVibeSeriesSuffixNumber = 14;
            variableValues{iVariable} = {sprintf('RADIAL_VIBE_SERIES_SUFFIX %025.10f', radialVibeSeriesSuffixNumber)};            
        case 'TriggerTime'
            t = sprintf('TriggerTime %030.15f', dicomHeader.TriggerTime);
            variableValues{iVariable} = {t};
        case 'SeriesNumber'
            variableValues{iVariable} = dicomHeader.SeriesNumber;
        case 'FlipAngle'
            t = sprintf('FlipAngle %030.15f', dicomHeader.FlipAngle);
            variableValues{iVariable} = {t};
        case 'EchoTime'
            variableValues{iVariable} = dicomHeader.EchoTime;
        otherwise
            error('Uknown indexing variable "%s".', variableNames{iVariable});
    end
end

variableValues = table(variableValues{:},'VariableNames', variableNames);