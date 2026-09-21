function [contentLoaded, content, contentTimestamp] = loadIfExistAndTimestampNewerThan(contentFile, referenceTimestamp)

contentLoaded = false;
content = [];
contentTimestamp = [];

if exist(contentFile, 'file') == 2
    
    fileStructure = matfile(contentFile);
    
    contentTimestamp = fileStructure.timestamp;
    
    if iso8601ToDatetime(referenceTimestamp) < iso8601ToDatetime(contentTimestamp)
        
        content = fileStructure.content;
        contentLoaded = true;
        
    end
    
end