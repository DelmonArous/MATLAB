function [img, bw, tform] = getCellContainerBorder(path, img, fig)

% The following variables are required for execution:
%   path = directory path to where the cell container image is saved,
%       showing an image of an empty cell container identical to the
%       container used for cell colony cultivation.
%   img = single channel image matrix of the bordered cell flask/dish cell
%       inhabiting cell colonies.
%
% The following variables are returned upon successful completion when
% input arguments are provided:
%   bw = binary mask of the cell container borders.

%% Get filename
[~, filename , ~] = fileparts(path);

%% Read cell container image

% FOR BAKTERIER!!
% img = load(path);
% % % %     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Colony Counter paper\Results ICPR2020\Dish Template Masks\E coli masks\Ecoli_control2_mask.txt');
% % % % img = load( ...
% % % %     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\BW.txt');
% bw = double(~img); % double(~img) for T47D
% bw = (bw - min(min(bw))) ./ abs(max(max(bw)) - min(min(bw)));
% bw = imerode(bw, strel('disk', 4)); % strel('disk', 4)) for T47D
% % % % figure(); imshow(bw, [])
%
% [img_ref, ~] = imread(path_ref);
% img_ref = double(rgb2gray(img_ref));
%
% [optimizer, metric] = imregconfig('multimodal');
%
% tform = imregtform(img_fixed, img_ref, 'rigid', optimizer, metric);
% img_fixed = imwarp(img_fixed, tform, 'OutputView', imref2d(size(img_fixed)));
%
% tform = imregtform(img_pca, img_ref, 'rigid', optimizer, metric);
% img_pca = imwarp(img_pca, tform, 'OutputView', imref2d(size(img_fixed)));
%
% tform = imregtform(bw, img_fixed, 'rigid', optimizer, metric);
% bw = imwarp(bw, tform, 'OutputView', imref2d(size(img_fixed)));


% Attempt to load image file using imread
try
    
    % If imread is successful, store the image information
    [img_fixed_rgb, ~] = imread(path);
    
catch
    
    % Otherwise, the file is either corrupt or not a real image.
    % Throw an error and stop execution
    uialert(fig, ...
        ['File ', filename, ' is not a valid cell container image!'], ...
        'Invalid image file');
    return;
    
end

%% Extraction of image plane

% Dimensions of container image data array
[~, ~, dims] = size(img_fixed_rgb);

% Check if the image array contains RGB channels
if dims == 3
    
    % Convert the RGB (truecolor) image to grayscale
    img_fixed = double(rgb2gray(img_fixed_rgb));
    
elseif dims < 3
    
    % Otherwise, use the first (red) channel as the reference image
    img_fixed = double(img_fixed_rgb(:,:,1));
    
else
    
    % Throw an error and stop execution
    uialert(fig, ...
        'The container image data array contains too many dimensions.', ...
        'Invalid image array dimensions');
    return;
    
end

%% Min-max normalization between 0-1

% img_moving_gray = double(rgb2gray(img_moving));
% % img_moving_hsv  = double(rgb2hsv(img_moving));
% % img_moving_hsv  = img_moving_hsv(:,:,3);
% 
% % Normalization
% [img_moving_norm, ~]    = normalizationMinMax(img_moving_gray);
% [img_fixed_norm, ~]     = normalizationMinMax(img_fixed);

% img_moving_norm = img;
img_fixed = ( img_fixed - min(img_fixed(:)) ) ./ ...
    abs( max(img_fixed(:)) - min(img_fixed(:)) );

%% Estimate geometric transformation that aligns the two 2D images

% getCellContainerBorder(path_fixed, path, img_fixed, fig)
% % if ~strcmp(path_fixed, path)

[optimizer, metric] = imregconfig('multimodal');
tform = imregtform(img, img_fixed, 'rigid', optimizer, metric);
Rfixed = imref2d(size(img_fixed));
img = imwarp(img, tform, 'OutputView', Rfixed);

% tform = imregtform(img_fixed, img, 'rigid', optimizer, metric);
% Rfixed = imref2d(size(img_fixed));
% img_fixed = imwarp(img_fixed, tform, 'OutputView', Rfixed);



% Start try-catch block to safely test for CUDA functionality
%     try
%
%         % Clear and initialize GPU memory. If CUDA is not enabled, or if the
%         % Parallel Computing Toolbox is not installed, this will error, and the
%         % function will automatically rever to CPU computation via the catch
%         % statement
%         gpuDevice(1);
%
%         tform = gather(imregtform(gpuArray(bw), gpuArray(img_moving), ...
%             'rigid', optimizer, metric));
%         img_fixed_norm = gather(imwarp(gpuArray(bw), tform, 'OutputView', ...
%             imref2d(gpuArray(size(img_moving)))));
%
%     catch

% If GPU fails, revert to CPU computation
% tform = imregtform(img_moving_norm, img_fixed_norm, 'rigid', ...
%     optimizer, metric);
% Rfixed = imref2d(size(img_fixed_norm));
%         img_moving_norm = imwarp(img_moving_norm, tform, 'OutputView', ...
%             Rfixed);

%     end

% end

%% Cell dish/flask extraction

% Pre-processing of the image; Gaussian filtering and adaptive histogram
% equalization
img_fixed = imgaussfilt(img_fixed, [3 3]);
img_fixed = adapthisteq(img_fixed, 'NumTiles', [16 16], ...
    'ClipLimit', 0.02, 'Distribution', 'exponential');

if ~isempty(size(img_fixed, 3)) && size(img_fixed,3) == 3
    img_fixed = double(mean(img_fixed,3));
elseif ~isempty(size(img_fixed,3)) && size(img_fixed,3) == 2
    img_fixed = double(mean(img_fixed,2));
else
    img_fixed = double(img_fixed);
end

% Threshold detection using most occuring value in histogram
% [counts, x] = hist(img_fixed(:), 255);
% delta = 0.45; % 0.25
% bw1 = img_fixed <= x(find(counts == max(max(counts)))) + delta;
% bw2 = img_fixed >= x(find(counts == max(max(counts)))) - delta;
% bw = bw1 & bw2;

% Otsu's method to binarize template image
bw = imbinarize(img_fixed, graythresh(img_fixed));

figure(); imshow(bw, [])

%% Clean up by applying morphological operations
% bw = imclearborder(bw, 8);
bw(1:15,:)          = 0;
bw(end-15:end,:)    = 0;
bw = imdilate(bw, strel('disk', 4));
bw = bwmorph(bw, 'majority');
bw = bwmorph(bw, 'thicken');
bw = bwmorph(bw, 'bridge');

% figure(); imshow(bw, [])

bw = imerode(bw, strel('disk', 6));
bw = imclearborder(bw, 8);
bw = bwareaopen(bw, 20000);
bw = imfill(bw, 'holes');

figure(); imshow(bw, [])

bw = bwmorph(bw, 'clean');
bw = imopen(bw, strel('disk', 5));
bw = imclose(bw, strel('disk', 100)); % 100 eller 500 egentlig; 200!!! 40?
bw = imfill(bw, 'holes');
bw = bwareaopen(bw, 250000);

figure(); imshow(bw, [])

bw = imdilate(bw, strel('disk', 25)); % 12 % 15?
bw = bwmorph(bw, 'majority');
bw = bwmorph(bw, 'thicken');
bw(1:15,:)          = 0; 
bw(:,1:15)          = 0;
bw(end-15:end,:)    = 0; 
bw(:,end-15:end)    = 0;

figure(); imshow(bw, [])
figure(); imshow(img, [])

%% Check if the extracted region is too small
% while nnz(bw) < 0.65*size(bw,1)*size(bw,2) && delta > 0 % 0.65
%     
%     % If so, iteratively decrease the span around the most occuring value
%     delta = delta - 0.01; % 0.025
%     
%     % Threshold detection using most occuring value in histogram
%     bw1 = img_fixed <= x(find(counts == max(max(counts)))) + delta;
%     bw2 = img_fixed >= x(find(counts == max(max(counts)))) - delta;
%     bw = bw1 & bw2;
%     
%     % Clean up by applying morphological operations
%     bw(1:15,:) = 0;
%     bw(end-15:end,:) = 0;
%     bw = imerode(bw, strel('disk', 10));
%     bw = imclearborder(bw, 8);
%     bw = imdilate(bw, ones(3));
%     bw = bwmorph(bw, 'bridge', Inf);
%     bw = bwareaopen(bw, 16000);
%     bw = imfill(bw, 'holes');
%     bw = bwmorph(bw, 'clean', Inf);
%     bw = imopen(bw, strel('disk', 4));
%     bw = imclose(bw, strel('disk', 40));
%     bw = imfill(bw, 'holes');
%     bw = bwareaopen(bw, 250000);
% %     bw = imdilate(bw, ones(15));
%             
% end

%% Clear temporary variables
clear filename dims optimizer metric counts x bw1 bw2 delta;

end