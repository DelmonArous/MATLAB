function varargout = runPCA(varargin)

% The following variables are required for execution:
%   varargin{1} = multivariate image data composed of color channel(s) of a cell
%       dish containing stained colonies. This image will be reshaped into
%       a 2-D matrix X of size rows * cols-by-dims.
%   varargin{2} = principal component coefficients (EIGVEC) for control
%       colony image dishes in order transform treatment dishes (i.e.
%       tranfer learning).
%
% The following variables are returned upon successful completion when
% input arguments are provided:
%   varargout{1} = principal component coefficients (EIGVEC) for the
%       rows * cols-by-dims data matrix X.
%   varargout{2} = principal component score (SCORE), which represents
%       the transform into PCA space.
%   varargout{3} = principal component variances (LATENT), i.e., the
%       eigenvalues of the covariance matrix of X.
%   varargout{4} = the 1st principal component image.
%   varargout{5} = the 2nd principal component image.
%   varargout{6} = the 3rd principal component image.

%% Size and number of color channels
[rows, cols, dims] = size(varargin{1});

%% Get an N by dims array of all color plane values.

% Reshape image into rows*cols by dims 2-D array, so that each column
% contains pixel values from every image (color) channel
listOfChannelValues = double(reshape(varargin{1}, rows * cols, dims));

%% Principal Component Analysis (PCA)

if nargin == 2
    
    % Remove the mean variable-wise (row-wise)
    X = listOfChannelValues.';
    B = X - repmat(mean(X,2), 1, size(X,2));
    
%     % Calculate eigenvectors and eigenvalues of the covariance matrix
%     [eigvec, ~] = eig(cov(B'));
%     
%     % Order by largest eigenvalue
%     eigvec = eigvec(:, end:-1:1); eigvec = eigvec';
    eigvec = varargin{2}; eigvec = eigvec';

    % Generate PCA component space (PCA scores)
    score = eigvec * B; score = score'; % score(:,2) = (-1) .* score(:,2);
    
    % plot PCA space of the first two PCs: PC1 and PC2
    % figure();
    % plot(score(:,1),score(:,2), '.')
    
else
    
    [eigvec, score] = pca(listOfChannelValues);
    
end


if dims == 3
    
    % Column 1, 2 and 3 contains the values of 1st, 2nd and 3rd principal
    % component. Extract each column and reshape back into a rectangular
    % image of the same size as the input image
    img_pca1 = double(reshape(score(:,1), rows, cols));
    img_pca2 = double(reshape(score(:,2), rows, cols));
    img_pca3 = double(reshape(score(:,3), rows, cols));
    
else
    
    % Column 1 contains the values of 1st principal component. Extract
    % column and reshape back into a rectangular image of the same size
    % as the input image
    img_pca1 = double(reshape(score(:,1), rows, cols));
    
end

%% Store return variables

% Check if input image consists of RGB channels
if dims == 3
    
    if nargout >= 1; varargout{1} = eigvec; end
    if nargout >= 2; varargout{2} = score; end
    if nargout >= 3; varargout{3} = img_pca1; end
    if nargout >= 4; varargout{4} = img_pca2; end
    if nargout >= 5; varargout{5} = img_pca3; end
    
else
    
    if nargout >= 1; varargout{1} = eigvec; end
    if nargout >= 2; varargout{2} = score; end
    if nargout >= 3; varargout{3} = img_pca1; end
    
end

%% Clear temporary variables
clear rows cols dims listOfChannelValues X B eigvec score latent ...
    img_pca1 img_pca2 img_pca3;

end