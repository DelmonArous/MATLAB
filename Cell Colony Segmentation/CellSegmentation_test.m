function varargout = CellSegmentation_test(varargin)

% The following variables are required for execution:
%   varargin{1} = directory path to where the cell colony image is saved,
%       showing cell colonies inhabiting a bordered cell container
%       (flask/dish).
%   varargin{2} = directory path to where the container image is saved,
%       showing an image of an empty cell container identical to the
%       container used for cell colony cultivation. This image should
%       preferrably be of the same dimensions as the colony image. Also,
%       the container should be pictured in the same position as the
%       cell colony image to extract out the inner region of the container.
%   varargin{3} = directory path to where the results (plots, image arrays)
%       should be saved.
%   varargin{4} = parameters needed for segmentation
%   varargin{5} = user interface object inherited from App Designer
%
% The following variables are returned upon successful completion when
% input arguments are provided:
%   varargout{1} = cell staining density estimate of the inner cell
%       container region.
%   varargout{2} = property structure of cell colony area estimate of each
%       segmented colony in the image.
%   varargout{3} = property structure of median red pixel estimate of each
%       segmented colony in the image.
%   varargout{4} = property structure of median green pixel estimate of each
%       segmented colony in the image.
%   varargout{5} = property structure of median blue pixel estimate of each
%       segmented colony in the image.
%   varargout{6} = property structure of median gray pixel estimate of each
%       segmented colony in the image.
%   varargout{7} = property structure of median PCA1 pixel estimate of each
%       segmented colony in the image.
%   varargout{8} = property structure of median PCA2 pixel estimate of each
%       segmented colony in the image.
%   varargout{9} = property structure of median PCA3 pixel estimate of each
%       segmented colony in the image.

%% Verify at least four input arguments are provided
if nargin < 5
    
    % If not, throw an error and stop execution
    uialert(varargin{5}, ...
        'Too few input arguments passed into function.', ...
        'Invalid input arguments');
    return;
    
end

% Get filename and file extension
[~, filename , ext] = fileparts(varargin{1});

%% If a valid screen size is returned (MATLAB was run without -nodisplay)

if usejava('jvm') && feature('ShowFigureWindows')
    
    % Start progress bar
    progress = uiprogressdlg(varargin{5}, 'Title', ...
        'Segmenting colonies...', 'Message', '1', 'Cancelable', 'on');
    progress.Value = 0.0;
    progress.Message = sprintf( ...
        '%s%s:\nReading colony image', ...
        filename, ext);
    
end

%% Read cell colony image

% Attempt to load image file using imread
try
    
    % If imread is successful, store the image information
    [img_in, ~] = imread(varargin{1});
    img_original = img_in;
    
    % Store segmentation parameters
    param = varargin{4};
    
catch
    
    % Otherwise, the file is either corrupt or not a real image.
    % Throw an error and stop execution
    uialert(varargin{5}, ...
        ['File ', filename, ' is not a valid cell colony image!'], ...
        'Invalid image file');
    return;
    
end

%% Extraction of image planes

% Dimensions of input image data array
[rows, cols, dims] = size(img_in);

% Check if the image array contains RGB channels. If so, extract RGB color
% channels and normalize them between [0,1]
if dims == 3
    
    % Store red, green and blue channel, respectively
    redchannel = double(img_in(:,:,1));
    greenchannel = double(img_in(:,:,2));
    bluechannel = double(img_in(:,:,3));
        
    % Convert the RGB (truecolor) image to grayscale
    img_in = double(rgb2gray(img_in));
    
    % Min-max normalization of each color channel image
    img_red = (redchannel - min(min(redchannel))) ./ ...
        abs(max(max(redchannel)) - min(min(redchannel)));
    img_green = (greenchannel - min(min(greenchannel))) ./ ...
        abs(max(max(greenchannel)) - min(min(greenchannel)));
    img_blue = (bluechannel - min(min(bluechannel))) ./ ...
        abs(max(max(bluechannel)) - min(min(bluechannel)));
    
elseif dims < 3
    
    % Otherwise, use the first (red) channel as the reference image
    img_in = double(img_in(:,:,1));
    
else
    
    % Throw an error and stop execution
    uialert(varargin{5}, ...
        'The image data array contains too many dimensions.', ...
        'Invalid image dimensions');
    return;
    
end

%% Update progress bar

if exist('progress', 'var')
    
    % Check for Cancel button press
    if progress.CancelRequested
        close(progress);
        delete(progress);
        delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
        clear all;
        close all;
        return;
    end
    
    % Update progress
    progress.Value = 0.1;
    progress.Message = sprintf( ... 
        '%s%s:\nColor channel extraction and normalization', ...
        filename, ext);
    
end

%% Min-max normalization between 0-1

% Background check
if median(median(img_in)) < max(max(img_in))
    npx_white = numel(find(img_in >= median(median(img_in))));
    npx_black = numel(find(img_in < median(median(img_in))));
else
    npx_white = numel(find(img_in >= mean(mean(img_in))));
    npx_black = numel(find(img_in < mean(mean(img_in))));
end

% Normalization
img_norm2 = (img_in - min(min(img_in))) ./ ...
    abs(max(max(img_in)) - min(min(img_in)));
%img_norm2 = MinMaxNorm(img_in);

if ~isempty(npx_white) && ~isempty(npx_black) && (npx_white > npx_black) % mean(mean(img_norm2)) > 0.45
    img_in                      = imcomplement(double(img_in(:,:,1)));
    img_norm                    = (img_in - min(min(img_in))) / ...
        abs(max(max(img_in)) - min(min(img_in))); % MinMaxNorm(img_in);
    img_norm(img_norm > 0.9)    = 0.9;
else
    img_norm                    = img_norm2;
end

% [img_norm, img_norm2] = normalizationMinMax(img_in);

% Store normalized image (either grayscale or red channel)
img_plane = img_norm;

%% Update progress bar

if exist('progress', 'var')
    
    % Check for Cancel button press
    if progress.CancelRequested
        close(progress);
        delete(progress);
        delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
        clear all;
        close all;
        return;
    end
    
    % Update progress
    progress.Value = 0.2;
    progress.Message = sprintf( ...
        '%s%s:\nCell dish border extraction', ...
        filename, ext);
    
end

%% Attain only the inner area of the cell dish/flask

% Mask out the inner space within the dish/flask
bw_border = getCellContainerBorder(varargin{1}, varargin{2}, ...
    img_norm, varargin{5});

% Multiply the dish/flask mask with the colony image to extract out
% solely free viable area in the container
img_norm = bw_border .* img_norm;

%% Update progress bar

if exist('progress', 'var')
    
    % Check for Cancel button press
    if progress.CancelRequested
        close(progress);
        delete(progress);
        delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
        clear all;
        close all;
        return;
    end
    
    % Update progress
    progress.Value = 0.3;
    progress.Message = sprintf( ...
        '%s%s:\nPrincipal component analysis', ...
        filename, ext);
    
end

%% Principal component analysis (PCA)

if dims == 3
    
    % Perform PCA on the three bundled color channels
    [~, ~, latent, img_pca1, img_pca2, img_pca3] = ...
        runPCA(double(img_original));
        
    latent(2)/sum(latent)
    
    % The clonogenic information is usually projected onto the
    % 2nd principal component plane
    if ~isempty(img_pca2) && (latent(2)/sum(latent)) > 0.01 % egentlig > 0.01
        
        % Therefore, store the PCA2 image as the primary PCA image
        img_pca = img_pca2;
        param.pca = 'PCA2';
        
    elseif (isempty(img_pca2) || (latent(2)/sum(latent)) < 0.01) && ...
            ~isempty(img_pca1)
        
        % Otherwise, not enough color variation is provided in the
        % input image such that no information is projected onto the
        % 2nd principal component plane, but rather onto the 1st
        % prinicpal component plane. Therefore, store the PCA1 image
        % as the PCA image
        img_pca = img_pca1;
        param.pca = 'PCA1';
        
    else
        
        % Throw an error and stop execution
        uialert(varargin{5}, ...
            'Principal component analysis (PCA) failed!', ...
            'Error in PCA acquisition');
        return;
                
    end
    
    % Otherwise, the input image consist of one
else
    
    % Perform PCA on the single image channel
    [~, ~, ~, img_pca1] = runPCA(double(img_original));
    
    % The clonogenic information is now projected onto the
    % 1st principal component plane
    if ~isempty(img_pca1)
        
        % Therefore, store the PCA1 image as the primary PCA image
        img_pca = img_pca1;
        param.pca = 'PCA1';
        
    else
        
        % Throw an error and stop execution
        uialert(varargin{5}, ...
            'Principal component analysis (PCA) failed!', ...
            'Error in PCA acquisition');
        return;
                
    end
    
end

%% Update progress bar

if exist('progress', 'var')
    
    % Check for Cancel button press
    if progress.CancelRequested
        close(progress);
        delete(progress);
        delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
        clear all;
        close all;
        return;
    end
    
    % Update progress
    progress.Value = 0.35;
    progress.Message = sprintf( ...
        '%s%s:\nFiltration and histogram equalization', ...
        filename, ext);
    
end

%% Guassian smoothing and removal of outliers

% Gaussian filtering
img_norm1 = imgaussfilt(img_norm, param.gaussfilt_size);

% Opening-by-Reconstruction
img_e = imerode(img_norm1, strel('disk', param.obrcbr_size));
img_obr = imreconstruct(img_e, img_norm1);

% Opening-Closing-by-Reconstruction
img_obrd = imdilate(img_obr, strel('disk', param.obrcbr_size));
img_obrcbr = imreconstruct(imcomplement(img_obrd), imcomplement(img_obr));
img_norm1 = imcomplement(img_obrcbr);
img_norm1(img_norm1 > 0.9) = 0.9;

%% Background correction by top-hat filtering

% img_norm = imsubtract(imadd(img_norm1,  ...
%     imtophat(img_norm1, strel('disk', 90))), ...
%     imbothat(img_norm1, strel('disk', 90)));
% [img_norm, ~] = normalizationMinMax(img_norm);

if ~isempty(npx_white) && ~isempty(npx_black) && (npx_white > npx_black) % mean(mean(img_norm2)) > 0.45
    img_norm = imbothat(img_norm1, strel('disk', param.tophat_size));
else
    img_norm = imtophat(img_norm1, strel('disk', param.tophat_size));
end

%% Perform CLAHE (histogram equalization)

maxthresh                           = prctile(double(img_norm(:)), 99);
minthresh                           = prctile(double(img_norm(:)), 1);
img_norm(img_norm > maxthresh)      = maxthresh;
img_norm(img_norm < minthresh)      = 0;
img_norm = (adapthisteq(img_norm) - min(min(adapthisteq(img_norm)))) ./ ...
    abs(max(max(adapthisteq(img_norm))) - min(min(adapthisteq(img_norm))));

figure(); imshow(img_norm, [])

%% Update progress bar

if exist('progress', 'var')
    
    % Check for Cancel button press
    if progress.CancelRequested
        close(progress);
        delete(progress);
        delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
        clear all;
        close all;
        return;
    end
    
    % Update progress
    progress.Value = 0.4;
    progress.Message = sprintf( ...
        '%s%s:\nBLOB extraction by K-means', ...
        filename, ext);
    
end

%% Mask out Binary Large OBjects (BLOBs) in input image by k-means

% Based on the defined search space, extract BLOBs of interest by creating
% a binary mask
bw_BLOBs = extractBLOBs(img_pca, bw_border, param);

%% Estimate staining density

area_total = sum(bw_border(:));
area_staining = sum(bw_BLOBs(:));

% Estimate cell density
confluency = area_staining/area_total;

%% Check if no BLOBs are identified

% If the BLOB mask contains only zero component values, then estimate
% BLOB features to return empty structures
if nnz(bw_BLOBs) == 0
    
    % Check if input image consists of RGB channels
    if dims == 3
        
        % If so, measure area along with RGB, grayscale, PCA1, PCA2 and
        % PCA3 distributions in the BLOB mask
        [areastats, redstats, greenstats, bluestats, graystats, ...
            pca1stats, pca2stats, pca3stats] = measureCellColonies( ...
            bw_BLOBs, img_red, img_green, img_blue, img_plane, ...
            img_pca1, img_pca2, img_pca3);
        
        % Store the return variables
        if nargout >= 1; varargout{1} = confluency; end
        if nargout >= 2; varargout{2} = areastats; end
        if nargout >= 3; varargout{3} = redstats; end
        if nargout >= 4; varargout{4} = greenstats; end
        if nargout >= 5; varargout{5} = bluestats; end
        if nargout >= 6; varargout{6} = graystats; end
        if nargout >= 7; varargout{7} = pca1stats; end
        if nargout >= 8; varargout{8} = pca2stats; end
        if nargout >= 9; varargout{9} = pca3stats; end
        
    else
        
        % Otherwise, a single image plane is provided. In this case,
        % measure area along with red and PCA1 distributions in the
        % BLOB mask
        [areastats, imgstats, pca1stats] = measureCellColonies( ...
            bw_BLOBs, img_plane, img_pca1);
        
        % Store the return variables
        if nargout >= 1; varargout{1} = confluency; end
        if nargout >= 2; varargout{2} = areastats; end
        if nargout >= 3; varargout{3} = imgstats; end
        if nargout >= 4; varargout{4} = pca1stats; end
        
    end
    
    % Update progress bar
    if exist('progress', 'var')
        
        % Update progress
        progress.Value = 1;
        progress.Message = 'Segmentation completed';
        
        % Check for Cancel button press
        if progress.CancelRequested
            close(progress);
            delete(progress);
            delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
            return;
        end
        
    end
    
    % Finish the image segmentation procedure
    return;
    
end

%% Update progress bar

if exist('progress', 'var')
    
    % Check for Cancel button press
    if progress.CancelRequested
        close(progress);
        delete(progress);
        delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
        clear all;
        close all;
        return;
    end
    
    % Update progress
    progress.Value = 0.5;
    progress.Message = sprintf( ...
        '%s%s:\nWatershed segmentation of BLOBs', ...
        filename, ext);
    
end

%% Segmentation of big BLOBs

% Define (global) watershed parameters
stats = regionprops(bw_BLOBs, img_norm, 'Area', 'Image');
param.medianarea = median([stats.Area]);

% Check if parallel computiation is neccessary, i.e. if there is a high
% BLOB density present that necessitates parallel acquisition for uniform
% distributed workload

if nnz(bw_BLOBs) > 0.01*rows*cols % 0.15
    
    nnz(bw_BLOBs) > 0.01*rows*cols
    
    % Try parallel acquisition to watershed segment the BLOBs
    try
        
        % Initialize parallel environment. If the Parallel Computing 
        % Toolbox is not installed/available, this will error, and the 
        % function will automatically revert to serial computation via 
        % the catch statement
        
        % If no pool, create new one
        poolobj = gcp;
        
        % Check if enough number of workers are available for a parallel
        % acquisition
        if poolobj.NumWorkers > 2
            
            % If so, split image into NumWorkers - 1 vertical stripes
            bw_BLOBs_par = [];
            for k = 1:(poolobj.NumWorkers-1)
                
                % Split input image into vertical stripes
                img_norm_temp{k} = img_norm(:, ...
                    floor( (k-1) * (cols/(poolobj.NumWorkers-1)) + 1): ...
                    floor( k * (cols/(poolobj.NumWorkers-1)) ));
                
                % Split corresponding BLOB-mask into vertical stripes
                bw_BLOBs_temp{k} = imclearborder(bw_BLOBs(:, ...
                    floor( (k-1) * (cols/(poolobj.NumWorkers-1)) + 1): ...
                    floor( k * (cols/(poolobj.NumWorkers-1)) )), 8);
                
                % Merge created binary stripes
                bw_BLOBs_par = horzcat(bw_BLOBs_par, bw_BLOBs_temp{k});
                
            end
            
            % Assign residual BLOBs to last worker
            img_norm_temp{poolobj.NumWorkers} = img_norm;
            bw_BLOBs_temp{poolobj.NumWorkers} = logical(bw_BLOBs - ...
                bw_BLOBs_par);
            
            % Parallel watershed acquisition of BLOBs
            parfor k = 1:poolobj.NumWorkers
                
                bw_seg_worker{k} = watershedBLOBs(img_norm_temp{k}, ...
                    bw_BLOBs_temp{k}, param, progress);
                
            end
            
            % Merge vertical watershed segmentation results from 
            % each worker
            bw_seg = [];
            for k = 1:poolobj.NumWorkers
                
                if k < poolobj.NumWorkers
                    bw_seg = logical(horzcat(bw_seg, bw_seg_worker{k}));
                elseif k == poolobj.NumWorkers
                    bw_seg = logical(bw_seg | bw_seg_worker{k});
                end
                
            end
            
        else
            
            % Serial acquisition if too few workers are available
            bw_seg = watershedBLOBs(img_norm, bw_BLOBs, param, progress);
            
        end
        
    catch
        
        % If parallell acquisition fails, revert to serial acquisition
        bw_seg = watershedBLOBs(img_norm, bw_BLOBs, param, progress);
        
    end
    
else
    
    % If there is a low density of BLOBs present, then perform serial 
    % acquisition for segmentation of BLOBs
    bw_seg = watershedBLOBs(img_norm, bw_BLOBs, param, progress);
    
end

%% Update progress bar

if exist('progress', 'var')
    
    % Check for Cancel button press
    if progress.CancelRequested
        close(progress);
        delete(progress);
        delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
        clear all;
        close all;
        return;
    end
    
    % Update progress
    progress.Value = 0.9;
    progress.Message = sprintf( ...
        '%s%s:\nPost-segmentation correction', ...
        filename, ext);
    
end

%% Post-segmentation correction

% Estimate mean pixel intensity of each colony
areastats = regionprops(bw_seg, img_norm, 'MeanIntensity');
all_meanintensities = [areastats.MeanIntensity];

% Remove segments that are too eccentric or have area and pixel
% intensity features that are out of defined range

bw_seg = bwpropfilt(bw_seg, 'Eccentricity', [0 0.98]); % 0.95 Hilde og 0.9 bakterier
bw_seg = bwpropfilt(bw_seg, 'Area', [param.area(1) param.area(4)]);

bw_seg = bwpropfilt(bw_seg, img_norm, 'MeanIntensity', ...
    [param.intensity_thresh*max(all_meanintensities(:)) ...
    max(all_meanintensities(:))]);

% Color-based segmentation using HSV space
% hsv = im2double(rgb2hsv(double(img_original)));
% hue = hsv(:,:,1) .* bw_seg;
% sat = hsv(:,:,2) .* bw_seg;
% temp_mask = hue >= 0.5 & hue <= 0.71; % aqua, teal, blue, navy, purple
% [r,c] = find(temp_mask);
% bw_seg = bwselect(bw_seg, c, r, 8);
% 
% satstats = regionprops(bw_seg, sat, 'MaxIntensity', 'PixelIdxList');
% for i = 1:numel(satstats)
%     sat(satstats(i).PixelIdxList) = 0.8*satstats(i).MaxIntensity;
% end
% 
% satstats = regionprops(bw_seg, sat, 'MeanIntensity');
% all_meansatintensities = [satstats.MeanIntensity];
% bw_seg = bwpropfilt(bw_seg, sat, 'MeanIntensity', ...
%     [param.intensity_thresh*max(all_meansatintensities(:)) ...
%     max(all_meansatintensities(:))]);

%% Colony feature estimations

% Check if input image consists of RGB channels
if dims == 3
    
    % If so, measure area along with RGB, grayscale, PCA1, PCA2 and PCA3
    % distributions in the segmented mask
    [areastats, redstats, greenstats, bluestats, graystats, ...
        pca1stats, pca2stats, pca3stats] = measureCellColonies( ...
        bw_seg, img_red, img_green, img_blue, img_plane, ...
        img_pca1, img_pca2, img_pca3);
    
else
    
    % Otherwise, a single image plane is provided. In this case, measure
    % area along with red and PCA1 distributions in the segmented mask
    [areastats, imgstats, pca1stats] = measureCellColonies( ...
        bw_seg, img_plane, img_pca1);
    
end

%% Plot delineated segmented colonies

h = figure();
% set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
imshow(img_original)

% Delineate segmented colonies
hold on
visboundaries(bw_seg, 'Color', 'red', 'LineWidth', 1)

% Add length scale bar
errorbar(0.1*cols, 0.95*rows, round(1/(param.px_size)), ...
    'horizontal', 'k.', 'LineWidth', 1.5, 'CapSize', 6);
text(0.1*cols, 0.95*rows, {'2 mm'}, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'Fontsize', 12);
hold off
set(gca, 'FontSize', 14)

% Save and close created plot
saveas(h, fullfile(varargin{3}, sprintf('%s-Seg.fig', filename)))
% saveas(h, fullfile(varargin{3}, sprintf('%s-Seg_%s_%s.png', ...
%     filename, num2str(param.gaussfilt_size(1)), num2str(param.obrcbr_size))))
% close(h)

%% Save final segmented binary colony mask
writematrix(bw_seg, fullfile(varargin{3}, ...
    sprintf('%s-SegMask.csv', filename)))
writematrix(bw_border, fullfile(varargin{3}, ...
    sprintf('%s-TemplateMask.csv', filename)))

%% Update progress bar

if exist('progress', 'var')
    
    % Check for Cancel button press
    if progress.CancelRequested
        close(progress);
        delete(progress);
        delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'));
        clear all;
        close all;
        return;
    end
    
    % Update progress
    progress.Value = 1.0;
    progress.Message = 'Segmentation completed';
    
end

%% Store return variables

% Check if input image consists of RGB channels
if dims == 3
    
    if nargout >= 1; varargout{1} = confluency; end
    if nargout >= 2; varargout{2} = areastats; end
    if nargout >= 3; varargout{3} = redstats; end
    if nargout >= 4; varargout{4} = greenstats; end
    if nargout >= 5; varargout{5} = bluestats; end
    if nargout >= 6; varargout{6} = graystats; end
    if nargout >= 7; varargout{7} = pca1stats; end
    if nargout >= 8; varargout{8} = pca2stats; end
    if nargout >= 9; varargout{9} = pca3stats; end
    
else
    
    if nargout >= 1; varargout{1} = confluency; end
    if nargout >= 2; varargout{2} = areastats; end
    if nargout >= 3; varargout{3} = imgstats; end
    if nargout >= 4; varargout{4} = pca1stats; end
    
end

%% Close and delete progress if it exists 

if exist('progress', 'var')
    
    close(progress);
    delete(progress);
    
end

end