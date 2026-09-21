function permissiveMakeDirIfNotExist(dirName, varargin)

if exist(dirName, 'dir') ~= 7
    
    permissiveMakeDir(dirName, varargin{:});
    
end

