function permissiveMakeDir(dirName, varargin)

groupName = [];

for iArgument = 1:2:numel(varargin)
    switch varargin{iArgument}
        case 'group'
            groupName = varargin{iArgument + 1};
        otherwise
            error('Unrecognized argument: %s', varargin{iArgument});
    end
end

commandLine = sprintf('mkdir -m a=rwx -p ''%s''', dirName);
[commandStatus, commandOutput] = system(commandLine);
assert(commandStatus == 0, '%s', commandOutput);

% if ~isempty(groupName)
%     commandLine = sprintf('chgrp ''%s'' ''%s''', groupName, dirName);    
%     [commandStatus, commandOutput] = system(commandLine);
%     assert(commandStatus == 0, '%s', commandOutput);
% end

