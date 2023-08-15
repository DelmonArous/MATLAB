function [img, bw_seg, ROI_img, ROI_seg, ROI_xy_GT, ROI_bb, N_seg, N_GT, ...
    precision, recall, F1, fn] = computeF1score(path_img, path_bw_seg, ...
    path_GT, dpi, ROI_size_mm)

% Read image
[~, fn, ~]  = fileparts(path_img);
img         = imread(path_img);

% Read segmentation mask and allocate respective centroids
bw_seg          = logical(readmatrix(path_bw_seg));
stats_seg       = regionprops(bw_seg, 'Centroid');
xycoord_seg     = round(cat(1, stats_seg.Centroid));

% Read ground truth (GT) xy-coordinates
T           = readtable(path_GT);
ROI_xy_GT   = round([T.XM T.YM]);
N_GT        = size(ROI_xy_GT, 1);

%     T               = readtable(filelist_test{i});
%     xycoord_segdata = [T.CentroidY_Coordinate_px_ T.CentroidX_Coordinate_px_];
%     xycoord_segdata = round(xycoord_segdata + 5.*randn(size(xycoord_segdata, 1), 2));

%% Compute ROI
px_size = 25.4/dpi;     % in mm/pixel
ROI_size_px = round(ROI_size_mm/px_size);

% Define scanning resolution, flask x- and y-coordinate range (in pixels)
x_range_px = [1 size(bw_seg, 1)];
y_range_px = [1 size(bw_seg, 2)];

% Define center of ROI
x_center = round((x_range_px(2) - x_range_px(1)) / 2);
y_center = round(((y_range_px(2) - y_range_px(1)) / 2) + y_range_px(1));

% Coordinate of upper left corner of ROI
x_bb = round(x_center - ROI_size_px(1)/2);
y_bb = round(y_center - ROI_size_px(2)/2);

% Define ROI bounding box
ROI_bb = [y_bb x_bb ROI_size_px(1)-1 ROI_size_px(2)-1];

% Crop and store the image segment which contains the selected
% region of interest (ROI)
ROI_seg = imcrop(bw_seg, ROI_bb);
ROI_seg = bwareaopen(ROI_seg, 15);
ROI_seg = imdilate(ROI_seg, strel('disk', 3));
ROI_seg = bwmorph(ROI_seg, 'thicken');
ROI_img = imcrop(img, ROI_bb);

% figure
% imshow(ROI_img)
% title(fn)
% hold on
% visboundaries(ROI_seg,'Color', 'black', 'LineWidth', 1)
% plot(ROI_xy_GT(:,1), ROI_xy_GT(:,2), 'b*')
% hold off

%% Find all segmented centroid coordinates located within ROI
%     xInRange    = (xycoord_seg(:,1) >= ROI_bb(1)) & ...
%         (xycoord_seg(:,1) <= ROI_bb(1) + ROI_bb(3));
%     yInRange    = (xycoord_seg(:,2) >= ROI_bb(2)) & ...
%         (xycoord_seg(:,2) <= ROI_bb(2) + ROI_bb(4));
%     xyInROI     = xInRange & yInRange;
%     xycoord_seg_cropped = xycoord_seg(xyInROI, :);
%     xycoord_seg_cropped(:,1) = xycoord_seg_cropped(:,1) - y_bb + 1;
%     xycoord_seg_cropped(:,2) = xycoord_seg_cropped(:,2) - x_bb + 1;

%% Assign a predicted colony number to each GT ROI number

% Label the predicted, binary colony image (ROI)
[L, N_seg] = bwlabel(ROI_seg);

% Assignment of GT coordinate positions
ind = sub2ind(size(L), ROI_xy_GT(:,2), ROI_xy_GT(:,1));
assign = L(ind);
assign = unique(assign(assign > 0)); % remove 0s (indicating that no predicted colony is at GT coordiantes) and duplicates

% Define numbers of TP, FN, FP
TP = length(assign);                        % number of uniquely assigned GT / predicted ROIs
FP = N_seg - length(assign);                % difference between total number of predicted ROIs and assigned GT
FN = size(ROI_xy_GT,1) - length(assign);    % difference between total number of GT ROIs and assigned GT

% Compute precision and recall
precision   = TP / (TP + FP);
recall      = TP / (TP + FN);

% Compute F1 score
if ~isnan(recall)
    F1      = 2 .* (precision * recall) ./ (precision + recall);
else
    recall  = 0;
    F1      = 0;
end

end