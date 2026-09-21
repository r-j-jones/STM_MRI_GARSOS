function timestamp = saveWithTimestamp(contentFile, content, timestamp) %#ok

%save(contentFile, '-v7.3', 'timestamp', 'content');

codeTimestamp = '2017-11-14 11:12:46 -05:00'; %#ok

lastwarn('');

save(contentFile, '-v7', 'timestamp', 'content', 'codeTimestamp');

[~, msgid] = lastwarn;

if strcmp(msgid, 'MATLAB:save:sizeTooBigForMATFile')
    
    warning('Failed to save content with v7 format. Trying with v7.3.');
    
    save(contentFile, '-v7.3', 'timestamp', 'content', 'codeTimestamp');
    
end





