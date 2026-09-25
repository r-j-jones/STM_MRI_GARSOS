function parms = setCustomSTMOptions(  )
% setDefaultSTMOptions manually sets the options for the STM algorithm.
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




%% Set the values manually here

%% STM computation parameters

tau      = 3;                         % Kernel radius. Default: 3

threshold = 0.035;                     % Threshold for C-matrix singular values. Default: 0.05
                                        % Note: In this example we don't use the default value.

M = 20;                               % Number of iterations for Orthogonal Iteration. Default: 30
                                        % Note: In this example we use a smaller value
                                        % to speed up the calculations.

interp_zp = 192;                       % Amount of zero-padding to create the low-resolution grid
                                        % if FFT-interpolation is used. Default: 24

gauss_win_param = 100;                % Parameter for the Gaussian apodizing window used to
                                        % generate the low-resolution image in the FFT-based
                                        % interpolation approach. This is the reciprocal of the
                                        % standard deviation of the Gaussian window. Default: 100

sketch_dim = 300;                     % Dimension of the sketch matrix used to calculate a
                                        % basis for the nullspace of the C matrix using a sketched SVD.
                                        % Default: 500. Note: In this example we use a smaller value
                                        % to speed up the calculations.

visualize_C_matrix_sv = 1;            % Binary variable. 1 = Singular values of the C matrix are displayed.
                                        % Default: 0.
                                        % Note: In this example we set it to 1 to visualize the singular values
                                        % of the C matrix. If sketched_SVD = 1 and if the curve of the singular values flattens out,
                                        % it suggests that the sketch dimension is appropriate for the data.
                                        
%% PISCO techniques

% The following techniques are used if the corresponding binary variable is equal to 1

kernel_shape = 1;                     % Binary variable. 1 = ellipsoidal shape is adopted for
                                        % the calculation of kernels (instead of rectangular shape).
                                        % Default: 1

FFT_nullspace_C_calculation = 1;      % Binary variable. 1 = FFT-based calculation of nullspace
                                        % vectors of C by calculating C'*C directly (instead of
                                        % calculating C first). Default: 1

sketched_SVD = 1;                     % Binary variable. 1 = sketched SVD is used to calculate
                                        % a basis for the nullspace of the C matrix (instead of
                                        % calculating the nullspace vectors directly and then the
                                        % basis). Default: 1

OrthogonalIteration_G_nullspace_vectors = 0; % Binary variable. 1 = Orthogonal Iteration approach is
                                        % used to find nullspace vectors of the G matrices
                                        % (instead of using SVD). Default: 1

FFT_interpolation = 1;                % Binary variable. 1 = sensitivity maps are calculated on
                                        % a small spatial grid and then interpolated to a grid with
                                        % nominal dimensions using an FFT-approach. Default: 1

verbose = 1;                          % Binary variable. 1 = PISCO information is displayed.
                                        % Default: 1


%% Now store the manual parameters in a structure
parms = struct();

parms.tau = tau;
parms.threshold = threshold;
parms.kernel_shape = kernel_shape;
parms.FFT_nullspace_C_calculation = FFT_nullspace_C_calculation;
parms.OrthogonalIteration_G_nullspace_vectors = OrthogonalIteration_G_nullspace_vectors;
parms.M = M;
parms.FFT_interpolation = FFT_interpolation;
parms.interp_zp = interp_zp;
parms.gauss_win_param = gauss_win_param;
parms.sketched_SVD = sketched_SVD;
parms.sketch_dim = sketch_dim;
parms.visualize_C_matrix_sv = visualize_C_matrix_sv;
parms.verbose = verbose;

end


