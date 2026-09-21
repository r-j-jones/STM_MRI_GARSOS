function exportArrayVolumeToNii(arrayVolume, dstFile)

fprintf('exportArrayVolumeToNii: %s\n', dstFile);

% niftiVolume = make_nii(arrayVolume.A, arrayVolume.v);

% A = [arrayVolume.R, arrayVolume.r0];
%
% niftiVolume.hdr.hist.sform_code = 1;
% niftiVolume.hdr.hist.srow_x = A(1,:);
% niftiVolume.hdr.hist.srow_y = A(2,:);
% niftiVolume.hdr.hist.srow_z = A(3,:);

% q = rotationmatrix_to_quaternion(arrayVolume.R);
% niftiVolume.hdr.hist.qform_code = 1;
% niftiVolume.hdr.hist.quatern_b = q(2);
% niftiVolume.hdr.hist.quatern_c = q(3);
% niftiVolume.hdr.hist.quatern_d = q(4);
% niftiVolume.hdr.hist.qoffset_x = arrayVolume.r0(1);
% niftiVolume.hdr.hist.qoffset_y = arrayVolume.r0(2);
% niftiVolume.hdr.hist.qoffset_z = arrayVolume.r0(3);

% save_nii(niftiVolume, dstFile);

A = [arrayVolume.R*diag(arrayVolume.v), arrayVolume.r0; 0 0 0 1];
B = A;
B(:,4) = -A*[1;1;1;1]+2*A(:,4);

% V.fname = dstFile;
nDimensions = max(ndims(arrayVolume.A),3);
arraySize = cell(1, nDimensions);
[arraySize{:}] = size(arrayVolume.A);
arraySize = cell2mat(arraySize);
% V.dim = arraySize(1:3);
% V.dt = [64 0];
% V.mat = B;
% V.pinfo = [1;0;0];
%
% V2 = spm_create_vol(V);
%
% spm_write_vol(V2, double(arrayVolume.A));

niftiVolume = nifti();
niftiVolume.dat = file_array(dstFile, arraySize, 'FLOAT32-LE');
niftiVolume.mat = B;
niftiVolume.mat0 = B;

% FIAT compatibility

niftiVolume.descrip = reshape(char(typecast(single(arrayVolume.r0),'uint8')),1,[]);
niftiVolume.aux_file = reshape(char(typecast(single(arrayVolume.R(1:6)),'uint8')),1,[]);

% Copy timestep if present

if iscell(arrayVolume.meta) && numel(arrayVolume.meta) == 1 
    
    t = 0;
    
    if isfield(arrayVolume.meta{1},'PixelDimensions') && numel(arrayVolume.meta{1}.PixelDimensions) >= 4
        t = arrayVolume.meta{1}.PixelDimensions(4);
    end
    
    if isfield(arrayVolume.meta{1},'pixdim') && numel(arrayVolume.meta{1}.pixdim) >= 5
        t = arrayVolume.meta{1}.pixdim(5);
    end
    
    niftiStructure = struct(niftiVolume);
    niftiStructure.hdr.pixdim(5) = t;
    niftiVolume = nifti(niftiStructure);
   
end

create(niftiVolume);
niftiVolume.dat(:) = arrayVolume.A(:);
















