function varargout = measureCellColonies(varargin)

% The following variables are required for execution:
%   varargin{1} = binary mask containing the segmented cell colonies.
%   varargin{>=2} = single channel image(s) (color, gray, PCA1, etc)
%       representing the cell colonies.
%
% The following variables are returned upon succesful completion when input
% arguments are provided:
%   varargout{1} = structure containing feature area measurements of the
%       segmented cell colonies provided in the binary mask varargin{1}.
%   varargout{>=2} = structure(s) containing feature measurements of the
%       segmented cell colonies provided in the binary mask varargin{1},
%       where the measurements (median and SD of pixel intensity) are based
%       upon the respectively given single channel image(s).

%% Measure area of image regions contained in the binary mask
areastats = regionprops(varargin{1}, 'Area', 'Centroid', 'Circularity', ...
    'Eccentricity', 'PixelIdxList');

% Store area measurements
varargout{1} = areastats;

%% Measure median pixel values and SD of segmented image regions

% Loop through each image plane provided
for i = 1:nargin-1
    
    % Extract and store pixel distributions of segmented regions in each
    % image plane
    varargout{i+1} = regionprops(varargin{1}, varargin{i+1}, 'PixelValues');
    
    % Check if there are any logical regions in the mask
    if ~isempty(areastats)
        
        % Loop through each binary segment (number of cell colonies)
        for j = 1:length(areastats)
            
            % Measure the median pixel value of each colony
            varargout{i+1}(j).MeanIntensity = ...
                mean(varargout{i+1}(j).PixelValues);
            
            % Measure corresponding SD of each binary region
            varargout{i+1}(j).SD = std(varargout{i+1}(j).PixelValues);
            
        end
        
    else
        
        % If not, return zero value in the first entry of the structure
        varargout{i+1}(1).MeanIntensity     = 0;
        varargout{i+1}(1).SD                = 0;
        
    end
    
end

end