function arrayVolume = loadArrayVolumeFromNii(srcFile)

% niftiVolume = load_nii(srcFile);
% 
% if niftiVolume.hdr.hist.sform_code ~= 0
% 
% A = [niftiVolume.hdr.hist.srow_x; niftiVolume.hdr.hist.srow_y; niftiVolume.hdr.hist.srow_z];
% 
% R = A(1:3,1:3);
% r0 = A(1:3,4);
% 
% else
%     R = eye(3,3);
%     r0 = zeros(3,1);
% end
% 
% v = niftiVolume.hdr.dime.pixdim(2:4)';

niftiVolume = nifti(srcFile);

[U,~,V] = svd(niftiVolume.mat0(1:3,1:3));
R = U*V';

v = sqrt(sum(niftiVolume.mat0(1:3,1:3).^2,1))';
r0 = niftiVolume.mat0*[1;1;1;1];
r0 = r0(1:3);

%A = spm_read_vols(spmVolume);

arrayVolume = ArrayVolume(double(niftiVolume.dat), R, v, r0, {niftiVolume.hdr});

