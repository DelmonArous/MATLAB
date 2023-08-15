function [I_rgb, I_red, I_green, I_blue, I_gray] = readEBT3film(path)

% Get filename and file extension
[~, filename , ~] = fileparts(path);

%% Read EBT3 film image

% Attempt to load image file using imread
try
    
    % If imread is successful, store the image information
    [I_rgb, ~] = imread(path);
    
catch
    
    % Otherwise, the file is either corrupt or not a real image.
    % Throw an error and stop execution
    error(['File ', filename, ' is not a valid cell colony image!']);
    
end

%% Extraction of image planes

% Dimensions of input image data array
[~, ~, dims] = size(I_rgb);

% Check if the image array contains RGB channels. If so, extract RGB color
% channels and normalize them between [0,1]
if dims == 3 || dims == 4
    
    % Store red, green and blue channel, respectively
    I_red = double(I_rgb(:,:,1));
    I_green = double(I_rgb(:,:,2));
    I_blue = double(I_rgb(:,:,3));
    I_rgb = I_rgb(:,:,1:3);
    
    % Convert the RGB (truecolor) image to grayscale
    I_gray = double(rgb2gray(I_rgb));
    
elseif dims < 3
    
    % Otherwise, use the first (red) channel as the reference image
    I_rgb = double(I_rgb(:,:,1));
    
else
    
    % Throw an error and stop execution
    error('The image data array contains too many dimensions.');
    
end

end
