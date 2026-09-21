% Class to represent medical 3d images

classdef ArrayVolume
    
    properties (Access  = public)
        A = [];
        R = [];
        v = []; % the size of a voxel, in mm
        r0 = [];
        meta = [];
    end    
    
    methods
        
        function obj = ArrayVolume(varargin)            
            if nargin >= 1; obj.A = varargin{1}; else obj.A = zeros(0, 0); end;
            if nargin >= 2; obj.R = varargin{2}; else obj.R = eye(3, 3); end;
            if nargin >= 3; obj.v = varargin{3}; else obj.v = ones(3, 1); end;
            if nargin >= 4; obj.r0 = varargin{4}; else obj.r0 = zeros(3, 1); end;
            if nargin >= 5; obj.meta = varargin{5}; else obj.meta = []; end;
        end       
        
        function nChannels = getNChannels(obj)
            arraySize = cell(1, 4);
            [arraySize{:}] = size(obj.A);
            nChannels = arraySize{4};
        end
        
        function arrayVolume = flatten(obj)
            
            arraySize = cell(1, 4);
            [arraySize{:}] = size(obj.A);
            arraySize = cell2mat(arraySize);
            
            if arraySize(4) == 1 && isreal(obj.A)
                arrayVolume = obj;
            else
                B = sqrt(sum(reshape(abs(obj.A), arraySize).^2,4));
                arrayVolume = ArrayVolume(B, obj.R, obj.v, obj.r0, obj.meta);
            end
            
        end
        
        function arrayVolumes = splitComplexComponents(obj)
            A1 = real(obj.A);
            A2 = imag(obj.A);
            arrayVolumes = cell(2, 1);
            arrayVolumes{1} = ArrayVolume(A1, obj.R, obj.v, obj.r0, obj.meta);
            arrayVolumes{2} = ArrayVolume(A2, obj.R, obj.v, obj.r0, obj.meta);
        end
        
        function arrayVolumes = splitChannels(obj)
            nDimensions = ndims(obj.A);
            nDimensions = max(nDimensions, 5);
            %arraySize = cell(1, nDimensions);
            %[arraySize{:}] = size(obj.A);
            %arraySize = cell2mat(arraySize);
            arraySize = size(obj.A);
            arraySize = num2cell(arraySize);
            arraySize(4:end) = cellfun(@(x)ones(1,x),arraySize(4:end),'UniformOutput',false);
            %channelSizes = arrayfun(@(x)ones(1,x),arraySize(4:end),'UniformOutput', false);
            B = mat2cell(obj.A, arraySize{:});
            B = permute(B, [4:nDimensions 1:3]);
            C = cellfun(@(x)ArrayVolume(x, obj.R, obj.v, obj.r0, obj.meta), B, 'UniformOutput',false);
            arrayVolumes = C; 
        end
        
        function b = cloneWithNewContent(obj, B)
            b = ArrayVolume(B, obj.R, obj.v, obj.r0, []);
        end
        
    end
    
    methods(Static)
        
        function arrayVolume = mergeChannels(arrayVolumes)
            
            B = cellfun(@(x)x.A, arrayVolumes, 'UniformOutput', false);
            
            n = max(ndims(arrayVolumes{1}.A), 3);
            
            m = ndims(arrayVolumes);
            
            B = permute(B, [(m + 1:m + n), 1:m]);
            
            B = cell2mat_mem(B);            
                        
            arrayVolume = ArrayVolume(B, arrayVolumes{1}.R, arrayVolumes{1}.v, arrayVolumes{1}.r0, arrayVolumes{1}.meta);
            
        end
        
        function arrayVolume = mergeChannelsDouble(arrayVolumes)
            for iVolume = 1:numel(arrayVolumes)
                arrayVolumes{iVolume}.A = double(arrayVolumes{iVolume}.A);
            end
            
            arrayVolume = ArrayVolume.mergeChannels(arrayVolumes);
        end
        
        function arrayVolume = mergeComplexComponents(arrayVolumes)
            
            assert(numel(arrayVolumes) == 2, 'Input must have two components');
            
            B = complex(arrayVolumes{1}.A, arrayVolumes{2}.A);
            
            arrayVolume = ArrayVolume(B, arrayVolumes{1}.R, arrayVolumes{1}.v, arrayVolumes{1}.r0, arrayVolumes{1}.meta);
            
        end
        
    end
    
end