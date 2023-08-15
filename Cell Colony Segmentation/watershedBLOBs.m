function bw_seg = watershedBLOBs(img, bw, param, progress)

% The following variables are required for execution:
%   img = input image that is either a grayscale or a single color (red)
%       channel image containing all inquired cell colony conglomerations
%       (BLOBs).
%   bw = binary mask containing the associated logical BLOBs for
%       watershed segmentation.
%   param = structure containing neccessary watershed parameters, such as
%       user-defined colony area, segmentation threshold and watershed
%       interval vector.
%
% The following variables are returned upon successful completion when
% input arguments are provided:
%   bw_seg = binary mask containing the segmented cell colonies.

%% Parameter extraction

% Feature properties of BLOBs
stats = regionprops(bw, img, 'Area', 'Image', 'PixelIdxList', ...
    'BoundingBox');

% Select only the BLOBs with sufficient area for watershed segmentation
param.thresh = 0.6*param.medianarea;            % lower area threshold
BLOBind = find([stats.Area] > param.thresh);    % BLOB indices
param.bb_buffer = 8;                            % boundingbox in pixels

% Initialization
watershed_temp = false(size(img));
bw_seg = logical(bw);

%% Watershed on each individual binary section of the mask

if (param.intensiveseg)
    
    for i = 1:length(BLOBind)
        
        % Update progress bar
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
            progress.Value = 0.5 + (0.89-0.5)*(i/length(BLOBind));
            
        end
        
        % Extract each BLOB by cropping (with a buffer 'safety' border) the
        % bounding box of each BLOB
        r = round(stats(BLOBind(i)).BoundingBox);
        bb_BLOB = [max(r(1)-param.bb_buffer,1) max(r(2)-param.bb_buffer,1) ...
            r(3)+2*param.bb_buffer-1 r(4)+2*param.bb_buffer-1];
        
        % Crop (with buffer border) the segment which contains the
        % selected BLOB
        img_BLOB = imcrop(img, bb_BLOB);
        
        % Crop the same segment accordingly from the binary mask
        bw_BLOB = false(size(bw));
        bw_BLOB(stats(BLOBind(i)).PixelIdxList) = 1;
        bw_BLOB = imcrop(bw_BLOB, bb_BLOB);
        
        % Watershed segmentation
        L = watershedBLOB(img_BLOB, bw_BLOB, param, param.ws_vec, 0);
        
        % Recursively watershed logical BLOB if segmented BLOB is
        % still too big
        temp_stats = regionprops(L, img_BLOB, 'Area');
        if ~isempty(find([temp_stats.Area] > param.thresh)) && ...
                max(max(bwlabel(L))) > max(max(bwlabel(bw_BLOB)))
            
            L = recursiveWatershedBLOB(img_BLOB, L, param);
            
        end
        
        %         img_original_cropped =  img_original( ...
        %             boundingboxBLOB(2):(boundingboxBLOB(2) + boundingboxBLOB(4)), ...
        %             boundingboxBLOB(1):(boundingboxBLOB(1) + boundingboxBLOB(3)), :);
        %         figure(); imshow(img_original_cropped)
        %         hold on
        %         visboundaries(L, 'Color', 'red', 'LineWidth', 1.5)
        %         hold off
        
        % Check if the segmented colony size exceeds the defined
        % area threshold
        temp_stats = regionprops(L, 'Area', 'Eccentricity', 'Circularity');
        if max(max(bwlabel(L))) == 1 && ...
                (temp_stats.Area > param.area(3)/2 || ...
                temp_stats.Circularity < 0.025)
            
            % If so, watershed segment the colony one final time
            L = watershedBLOB(img_BLOB, L, param, ...
                min(param.ws_vec)-0.03:0.01:0.45, 1);
            
        end
        
        % Replace inspected BLOB with its watershed segmented colonies
        if max(max(bwlabel(L))) > 0 && sum(sum(L)) >= 0.6*sum(sum(bw_BLOB))
            
            [f1, f2] = find(L > 0);
            fr = f1 + bb_BLOB(2) - 1;
            fc = f2 + bb_BLOB(1) - 1;
            watershed_temp(sub2ind(size(img), fr, fc)) = ...
                L(sub2ind(size(L), f1, f2));
            bw_seg(stats(BLOBind(i)).PixelIdxList) = 0;
            
        end
        
    end
    
else
    
    for i = 1:length(BLOBind)
        
        % Update progress bar
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
            progress.Value = 0.5 + (0.89-0.5)*(i/length(BLOBind));
            
        end
        
        % Extract each BLOB by cropping (with a buffer 'safety' border) the
        % bounding box of each BLOB
        r = round(stats(BLOBind(i)).BoundingBox);
        bb_BLOB = [max(r(1)-param.bb_buffer,1) max(r(2)-param.bb_buffer,1) ...
            r(3)+2*param.bb_buffer-1 r(4)+2*param.bb_buffer-1];
        
        % Crop (with buffer border) the segment which contains the
        % selected BLOB
        img_BLOB = imcrop(img, bb_BLOB);
        
        % Crop the same segment accordingly from the binary mask
        bw_BLOB = false(size(bw));
        bw_BLOB(stats(BLOBind(i)).PixelIdxList) = 1;
        bw_BLOB = imcrop(bw_BLOB, bb_BLOB);
        
        % Watershed segmentation
        L = watershedBLOB(img_BLOB, bw_BLOB, param, param.ws_vec, 0);
        
        % Replace inspected BLOB with its watershed segmented colonies
        if max(max(bwlabel(L))) > 0 && sum(sum(L)) >= 0.6*sum(sum(bw_BLOB))
            
            [f1, f2] = find(L > 0);
            fr = f1 + bb_BLOB(2) - 1;
            fc = f2 + bb_BLOB(1) - 1;
            watershed_temp(sub2ind(size(img), fr, fc)) = ...
                L(sub2ind(size(L), f1, f2));
            bw_seg(stats(BLOBind(i)).PixelIdxList) = 0;
            
        end
        
    end
    
end

bw_seg = imfill(logical(bw_seg) | watershed_temp, 'holes');
bw_seg = imerode(bw_seg, strel('disk', 3));

%% Clear temporary variables
clear stats param thresh BLOBind watershed_temp r ...
    bb_BLOB img_BLOB bw_BLOB L temp_stats f1 f2 fr fc

end