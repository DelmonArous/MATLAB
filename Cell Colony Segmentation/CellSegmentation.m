function varargout = CellSegmentation(varargin)

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
%   varargin{4} = parameters needed for segmentation.
%   varargin{5} = user interface object inherited from App Designer.
%
% The following variables are returned upon successful completion when
% input arguments are provided:
%   varargout{1} = cell staining density estimate of the inner cell
%       container region.
%   varargout{2} = property structure of cell colony area estimate of each
%       segmented colony in the image.
%   varargout{3} = property structure of median red pixel estimate of each
%       segmented colony in the image.
%   varargout{4} = property structure of median green pixel estimate of
%       each segmented colony in the image.
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
%   varargout{10} = principal component eigenvectors from the Principal
%       Component Analysis (PCA) needed for transfer learning.

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
    progress.Value = 0.05;
    progress.Message = sprintf( ...
        '%s%s:\nReading colony image', ...
        filename, ext);
    
end

%% Read cell colony image

% Attempt to load image file using imread
try
    
    % If imread is successful, store the image information
    [img_in, ~] = imread(varargin{1});
    img_original = img_in(:,:,1:3);
    
    % Dimensions of input image data array
    [rows, cols, dims] = size(img_in);
    
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
        '%s%s:\nPrincipal Component Analysis (PCA)', ...
        filename, ext);
    
end

%% Principal component analysis (PCA)

if dims == 3 || dims == 4
    
    % Perform PCA on the three bundled color channels
    [eigvec, Z_pca, pca1channel, pca2channel, pca3channel] = ...
        runPCA(double(img_original));

    % Throw an error and stop execution
    if isempty(pca1channel) && isempty(pca2channel)
        
        uialert(varargin{5}, ...
            'Principal component analysis (PCA) failed!', ...
            'Error in PCA acquisituon');
        return;
        
    end
    
    % Otherwise, the input image consist of one
else
    
    % Perform PCA on the single image channel
    [eigvec, Z_pca, pca1channel] = runPCA(double(img_original));
    
    % Throw an error and stop execution
    if isempty(pca1channel)
        
        uialert(varargin{5}, ...
            'Principal component analysis (PCA) failed!', ...
            'Error in PCA acquisituon');
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
    progress.Value = 0.2;
    progress.Message = sprintf( ...
        '%s%s:\nExtraction and normalization of image color/PC planes', ...
        filename, ext);
    
end

%% Extraction of image color and PC planes

% Check if the image array contains RGB channels. If so, extract RGB color
% channels and normalize them between [0,1]
if dims == 3 || dims == 4
    
    % Store red, green and blue channel, respectively
    redchannel      = double(img_in(:,:,1));
    greenchannel    = double(img_in(:,:,2));
    bluechannel     = double(img_in(:,:,3));
    graychannel     = double(rgb2gray(img_in(:,:,1:3))); 
    
    % Min-max normalization of each color channel image between 0-1
    img_red     = ( redchannel - min(redchannel(:)) ) ./ ...
        abs( max(redchannel(:)) - min(redchannel(:)) );
    img_green   = ( greenchannel - min(greenchannel(:)) ) ./ ...
        abs( max(greenchannel(:)) - min(greenchannel(:)) );
    img_blue    = ( bluechannel - min(bluechannel(:)) ) ./ ...
        abs( max(bluechannel(:)) - min(bluechannel(:)) );
    img_gray    = ( graychannel - min(graychannel(:)) ) ./ ...
        abs( max(graychannel(:)) - min(graychannel(:)) );
    img_pca1    = ( pca1channel - min(pca1channel(:)) ) ./ ...
        abs( max(pca1channel(:)) - min(pca1channel(:)) );
    img_pca2    = ( pca2channel - min(pca2channel(:)) ) ./ ...
        abs( max(pca2channel(:)) - min(pca2channel(:)) );
    img_pca3    = (pca3channel - min(pca3channel(:)) ) ./ ...
        abs( max(pca3channel(:)) - min(pca3channel(:)) );
    
elseif dims < 3
    
    % Store red channel
    redchannel = double(img_in(:,:,1));
    
    % Min-max normalization of the red color channel image
    img_red = ( redchannel - min(redchannel(:)) ) ./ ...
        abs( max(redchannel(:)) - min(redchannel(:)) );
    img_pca1    = ( pca1channel - min(pca1channel(:)) ) ./ ...
        abs( max(pca1channel(:)) - min(pca1channel(:)) );
    
else
    
    % Throw an error and stop execution
    uialert(varargin{5}, ...
        'The image data array contains too many dimensions.', ...
        'Invalid image dimensions');
    return;
    
end

% Store and use 1st principal component image as a reference intensity
% image for segmentation
img_in = double(img_gray); % double(img_pca1); % double(rgb2gray(img_in(:,:,1:3))); %  

%% Normalization and background check
img_norm2 = normalizationMinMax(img_in);

if mean(img_norm2(:)) > 0.45
    img_norm = normalizationMinMax(imcomplement(double(img_in(:,:,1))));
    img_norm(img_norm > 0.9) = 0.9;
else
    img_norm = img_norm2;
end

% [img_norm, img_norm2] = normalizationMinMax(img_in);
% img_plane = img_norm;
% figure(); imshow(img_norm, [])
% figure(); imshow(img_norm2, [])

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
        '%s%s:\nCell dish border extraction', ...
        filename, ext);
    
end

%% Attain only the inner area of the cell dish/flask

% Mask out the inner space within the dish/flask
if param.extractcontainer == 1
    [img_norm, bw_border, tform] = getCellContainerBorder(varargin{2}, ... % img_norm, bw_border, tform
        img_norm, varargin{5});
else
    bw_border = ones(rows, cols);
end

figure(); imshow(bw_border, [])
figure(); imshow(img_norm, [])

% Multiply the dish/flask mask with the colony image to extract out
% solely free viable area in the container
img_norm = bw_border .* img_norm;

figure(); imshow(img_norm, [])

%% PC channel selection

% Channel selection by GLCM contrast criterion
% PCstruct = selectPrincipalComponentChannel(Z_pca, rows, cols, dims, bw_border, tform);
PCstruct.channelopt.img = img_pca2;
PCstruct.channelopt.img = imwarp(PCstruct.channelopt.img, tform, ...
    'OutputView', imref2d(size(img_norm)));

% figure(); imshow(img_pca2, [])

%% Image registration of the input image to template image
for i = 1:dims
    img_original(:,:,i) = imwarp(img_original(:,:,i), tform, ...
        'OutputView', imref2d(size(img_norm)));
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
        '%s%s:\nGray-Level Co-Occurrence Matrix (GLCM) computation', ...
        filename, ext);
    
end

%% Guassian smoothing and removal of outliers

% Gaussian filtering
img_norm = imgaussfilt(img_norm, param.gaussfilt_size);

% img_e = imerode(img_norm, strel('disk', 6));
% img_obr = imreconstruct(img_e, img_norm);
% % Opening-Closing-by-Reconstruction
% img_obrd = imdilate(img_obr, strel('disk', 6));
% img_obrcbr = imreconstruct(imcomplement(img_obrd), imcomplement(img_obr));
% img_norm = imcomplement(img_obrcbr);

img_norm(img_norm > 0.9)    = 0.9;
img_norm(img_norm < 0.01)   = 0.0;

%% Background correction
% if mean(img_norm2(:)) > 0.45
%     img_norm = imbothat(img_norm, strel('disk', 90));
% else
%     img_norm = imtophat(img_norm, strel('disk', 90));
% end
% figure(); imshow(img_norm, [])
 
%% Perform CLAHE (histogram equalization)
img_norm = adapthisteq(img_norm);
img_norm = (img_norm - min(img_norm(:)) ) ./ ...
    abs( max(img_norm(:)) - min(img_norm(:)) );
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
bw_BLOBs = extractBLOBs(PCstruct.channelopt.img, bw_border, param); % param.area(2)

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
    if dims == 3 || dims == 4
        
        % If so, measure area along with RGB, grayscale, PCA1, PCA2 and
        % PCA3 distributions in the BLOB mask
        [areastats, redstats, greenstats, bluestats, graystats, ...
            pca1stats, pca2stats, pca3stats] = measureCellColonies( ...
            bw_BLOBs, img_red, img_green, img_blue, img_gray, ...
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
        [areastats, redstats, pca1stats] = measureCellColonies( ...
            bw_BLOBs, img_red, img_pca1);
        
        % Store the return variables
        if nargout >= 1; varargout{1} = confluency; end
        if nargout >= 2; varargout{2} = areastats; end
        if nargout >= 3; varargout{3} = redstats; end
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
if nnz(bw_BLOBs) > 0.15*rows*cols
    
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

% Post-processing of segmentation mask
bw_seg = bwmorph(bw_seg, 'hbreak');
bw_seg = bwmorph(bw_seg, 'spur');
bw_seg = bwmorph(bw_seg, 'majority');
bw_seg = bwmorph(bw_seg, 'thicken');

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

figure(); imshow(bw_BLOBs, [])
figure(); imshow(bw_seg, [])

% Estimate mean pixel intensity of each colony
areastats = regionprops(bw_seg, img_norm, 'MeanIntensity');
all_meanintensities = [areastats.MeanIntensity];

% Remove segments that are too eccentric or have area and pixel
% intensity features that are out of defined range
% bw_seg = bwpropfilt(bw_seg, 'Eccentricity', [0 0.99]); % 0.95 Hilde % 0.98 egentlig
bw_seg = bwpropfilt(bw_seg, 'Area', [param.area(1) param.area(4)*2]);
bw_seg = bwpropfilt(bw_seg, img_norm, 'MeanIntensity', ...
    [param.intensity_thresh*max(all_meanintensities(:)) ...
    max(all_meanintensities(:))]);

figure(); imshow(bw_seg, [])

% figure(); imshow(bw_seg, [])
% title('Watershed (processed)')
%
% % Color-based segmentation using HSV space
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
if dims == 3 || dims == 4
    
    % If so, measure area along with RGB, grayscale, PCA1, PCA2 and PCA3
    % distributions in the segmented mask
    [areastats, redstats, greenstats, bluestats, graystats, ...
        pca1stats, pca2stats, pca3stats] = measureCellColonies( ...
        bw_seg, img_red, img_green, img_blue, img_gray, ...
        img_pca1, img_pca2, img_pca3);
    
else
    
    % Otherwise, a single image plane is provided. In this case, measure
    % area along with red and PCA1 distributions in the segmented mask
    [areastats, redstats, pca1stats] = measureCellColonies( ...
        bw_seg, img_red, img_pca1);
    
end

%% Plot delineated segmented colonies

h = figure();
set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
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
% saveas(h, fullfile(varargin{3}, sprintf('%s-Seg.png', filename)))
close(h)

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
if dims == 3 || dims == 4
    
    if nargout >= 1;    varargout{1}  = confluency; end
    if nargout >= 2;    varargout{2}  = areastats; end
    if nargout >= 3;    varargout{3}  = redstats; end
    if nargout >= 4;    varargout{4}  = greenstats; end
    if nargout >= 5;    varargout{5}  = bluestats; end
    if nargout >= 6;    varargout{6}  = graystats; end
    if nargout >= 7;    varargout{7}  = pca1stats; end
    if nargout >= 8;    varargout{8}  = pca2stats; end
    if nargout >= 9;    varargout{9}  = pca3stats; end
    if nargout >= 10;   varargout{10} = eigvec; end
    
else
    
    if nargout >= 1;    varargout{1}  = confluency; end
    if nargout >= 2;    varargout{2}  = areastats; end
    if nargout >= 3;    varargout{3}  = redstats; end
    if nargout >= 4;    varargout{4}  = pca1stats; end
    
end

%% Close and delete progress if it exists

if exist('progress', 'var')
    
    close(progress);
    delete(progress);
    
end

end