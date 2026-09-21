function path2add(varargin)
    for i = 1 : nargin
        libPath = fullfile('./', varargin{i});
        if exist(libPath, 'dir') == 7
            addpath(libPath);
        else
            addpath(fullfile('/RadOnc-MRI1/Student_Folder/zeyiren/YZ_data_8_15/', varargin{i}));
        end
    end
end
