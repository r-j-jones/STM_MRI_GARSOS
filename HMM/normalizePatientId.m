function outputPatientId = normalizePatientId(inputPatientId)

if iscell(inputPatientId)
    
    outputPatientId = cellfun(@(x)normalizePatientId(x), inputPatientId, 'UniformOutput', false);
    
else
    
    if isnumeric(inputPatientId)
        
        outputPatientId = num2str(inputPatientId);
        
    else
        
        outputPatientId = strtrim(inputPatientId);
        
        if all(outputPatientId >= '0') && all(outputPatientId <= '9')
            outputPatientId = num2str(str2double(outputPatientId));
        end
        
    end
    
end