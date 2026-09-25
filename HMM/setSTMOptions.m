function pp = setSTMOptions( varargin )
% setDefaultSTMOptions sets the default options for the STM algorithm.
%
%   --tau:                             2D case: Parameter (in Nyquist units) that determines the size of 
%                                               the k-space kernel. For a rectangular kernel, the size is 
%                                               (2*tau+1) x (2*tau+1). For an ellipsoidal kernel, it is 
%                                               the radius of the associated neighborhood. Default: 3.
%
%   --threshold:                       Specifies how small a singular value needs to be (relative 
%                                      to the maximum singular value) before its associated 
%                                      singular vector is considered to be in the nullspace of 
%                                      the C-matrix. Default: 0.05.
%
%   --kernel_shape:                    Binary variable. 0 = rectangular kernel, 1 = ellipsoidal 
%                                      kernel. Default: 1.
%
%   --FFT_nullspace_C_calculation:     Binary variable. 0 = nullspace vectors 
%                                      of C are calculated from C'*C by calculating C first. 
%                                      1 = nullspace vectors of C are calculated from C'*C 
%                                      directly using an FFT-based approach. Default: 1.
%
%   --OrthogonalIteration_G_nullspace_vectors: Binary variable. 0 = nullspace 
%                                         vectors of the G matrices are calculated using SVD. 
%                                         1 = nullspace vectors of the G matrices are calculated 
%                                         using a Orthogonal Iteration approach. Default: 1.
%
%   --M:                               Number of iterations used in the Orthogonal Iteration approach 
%                                      to calculate the nullspace vectors of the G matrices. 
%                                      Default: 30.
%
%   --FFT_interpolation:               Binary variable. 0 = no interpolation. 1 = 
%                                      FFT-based interpolation is used. Default: 1.
%
%   --interp_zp:                       Amount of zero-padding to create the low-resolution grid 
%                                      if FFT-interpolation is used. 
%                                      2D case: The grid has dimensions 
%                                               (N1_cal + interp_zp) x (N2_cal + interp_zp) x Nc. 
%                                      Default: 24.
%
%   --gauss_win_param:                 Parameter for the Gaussian apodizing window used to 
%                                      generate the low-resolution image in the FFT-based 
%                                      interpolation approach. This is the reciprocal of the 
%                                      standard deviation of the Gaussian window. Default: 100.
%
%   --sketched_SVD:                    Binary variable. 1 = sketched SVD is used to calculate 
%                                      a basis for the nullspace of the C matrix. Default: 1.
%
%   --sketch_dim:                      Dimension of the sketch matrix used to calculate a basis 
%                                      for the nullspace of the C matrix using a sketched SVD. 
%                                      Only used if sketched_SVD is enabled. Default: 500.
%
%   --visualize_C_matrix_sv:           Binary variable. 1 = Singular values of the C matrix are displayed.
%                                      Default: 0. 
%                                      Note: If sketched_SVD = 1 and if the curve of the singular values flattens out,
%                                      it suggests that the sketch dimension is appropriate for the data.
%
%   --verbose:                         Binary variable. 1 = display PISCO information, including 
%                                      which techniques are employed and computation times for 
%                                      each step. Default: 1.




% Set default values for optional parameters
p = inputParser;

addParameter(p, 'tau', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'threshold', 0.05, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'kernel_shape', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'FFT_nullspace_C_calculation', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'OrthogonalIteration_G_nullspace_vectors', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'M', 30, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'FFT_interpolation', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'interp_zp', 24, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'gauss_win_param', 100, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'sketched_SVD', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'sketch_dim', 500, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'visualize_C_matrix_sv', 0, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'verbose', 1, @(x) isnumeric(x) && isscalar(x));

if isempty(varargin)
    parse(p);
else
    parse(p, varargin{:});
end

pp = p.Results;

end 




