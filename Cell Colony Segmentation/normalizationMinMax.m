function [img_norm] = normalizationMinMax(img)

% The following variables are required for execution: 
%   img = single channel image that needs to be normalized. 
%
% The following variables are returned upon successful completion when 
% input arguments are provided:
%   img_norm = min-max normalized image of img, adjusted for staining
%   contrast.
%   img_norm2 = min-max normalized image of img, unadjusted for staining
%   contrast.

%% Min-max normalization of img
img_norm2 = ( img - min(img(:)) ) ./ abs( max(img(:)) - min(img(:)) );

%% Background check
px_flag = false;

if median(img_norm2(:)) < max(img_norm2(:))
    npx_white = numel(find(img_norm2 >= median(img_norm2(:))));
    npx_black = numel(find(img_norm2 < median(img_norm2(:))));
else
    npx_white = numel(find(img_norm2 >= mean(img_norm2(:))));
    npx_black = numel(find(img_norm2 < mean(img_norm2(:))));
end

if (npx_white > npx_black) && ~isempty(npx_white) && ~isempty(npx_black)
    px_flag = true;
else
    px_flag = false;
end

%% Check if mean pixel intensity is more than 0.45
if px_flag %  && mean(mean(img_norm2)) > 0.45
    
    % If so, invert the image so the stained colonies have higher
    % pixel intensity (are brightest)
    img = imcomplement(double(img(:,:,1)));
    img_norm = ( img - min(img(:)) ) ./ abs( max(img(:)) - min(img(:)) );
    
else
    
    img_norm = img_norm2;
    
end

%% Clear temporary variables
clear img;

end