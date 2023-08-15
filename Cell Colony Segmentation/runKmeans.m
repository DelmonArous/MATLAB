function [bw] = runKmeans(img)

% The following variables are required for execution: 
%   img = input image containing the cell colony conglomerations (BLOBs). 
%
% The following variables are returned upon successful completion when 
% input arguments are provided:
%   bw = binary mask containing the BLOBs of the inquired cell colonies. 

%% Design the feature matrix Z
new_img = zeros(size(img)+2);
new_img(2:end-1,2:end-1) = double(img);

% Image boundary coordinates without first/last row/column
inner_coord = [2 2; size(new_img,1)-1 size(new_img,2)-1];

% 9x2 matrix with 1 row for the relative shift to reach neighbors
[d2, d1] = meshgrid(-1:1, -1:1);
displacement = [d1(:) d2(:)];

% Cell array to store all 9 shifted images
temp = {};

for i = 1:size(displacement,1)
    
   % X-indices of the submatrix when shifted by d(i,1)
   Xcoord = (inner_coord(1,1):inner_coord(2,1)) + displacement(i,1);
   
   % Y-indices of the submatrix when shifted by d(i,2)
   Ycoord = (inner_coord(1,2):inner_coord(2,2)) + displacement(i,2);
   
   % Image matrix resulting from shift by displacement(i,)
   temp{i} = reshape(new_img(Xcoord, Ycoord), 1, []);
   
end

% Column-wise bind all 9 shifted images (as vectors) from the array
Z_feat = vertcat(temp{:}).';

%% K-means acquisition

% Number of clusters that will definitely be present
k = 2; % background/foreground

% Define background-/foreground centroid-cluster vector
start_mat = [repelem(0.2,9); repelem(0.85,9)]; % repelem(1,9); repelem(0.2,9)
% opts = statset('Display','final');

% for i = 1:10

% The feature matrix Z is now size(img,1)*size(img,2)-by-9
[idx, C] = kmeans(Z_feat, k, 'Distance', 'sqeuclidean', ...
    'MaxIter', 5000, 'Start', start_mat); %  'Replicates', 1, 'Options', opts 
img_Kmeans = reshape(idx, size(img,1), size(img,2));

% Generate binary images for background and foreground pixels
bw_0 = (img_Kmeans == min(img_Kmeans(:)));
bw_1 = (img_Kmeans == max(img_Kmeans(:)));

% stats_0 = regionprops(bw_0, 'EulerNumber', 'Centroid', 'Circularity', ...
%     'ConvexArea', 'ConvexHull', 'ConvexImage', 'Extent');
% stats_1 = regionprops(bw_1, 'EulerNumber', 'Centroid', 'Circularity', ...
%     'ConvexArea', 'ConvexHull', 'ConvexImage', 'Extent');

% figure(); imshow(bw_0, [])
% figure(); imshow(bw_1, [])

% eva = evalclusters(Z_feat, idx, 'DaviesBouldin')

% figure;
% plot(Z_feat(idx==1,1),Z_feat(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(Z_feat(idx==2,1),Z_feat(idx==2,2),'b.','MarkerSize',12)
% plot(C(:,1), C(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3) 
% legend('Cluster 1', 'Cluster 2', 'Centroids', 'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off

% Get pixels that are labeled as colonies 
% Problem 1: if the image contains too many cells!!
% Problem 2: if the image contains too much noise!!

% [maxValue, indexOfMaxValue] = max(img_Kmean);

% Measure the solidity of all the blobs.

% figure(); imshow(img_Kmean, [])

% stats = regionprops(bw_border, img_Kmean, 'Area', 'Image', 'PixelIdxList', ...
%     'BoundingBox');
% r = round(stats.BoundingBox);
% bb = [max(r(1),1) max(r(2),1) r(3)-1 r(4)-1];
% % Crop (with buffer border) the segment which contains the
% % selected BLOB
% img_Kmean_crop = imcrop(img_Kmean, bb);
% 
% figure(); imshow(img_Kmean_crop, [])


% bw_border(bw_border == 0) = 3;
% temp_img_Kmean = img_Kmean .* bw_border;
% temp_img_Kmean(temp_img_Kmean >= 3) = 3;
% 
% figure(); imshow(temp_img_Kmean, [])
% 
% measurements = regionprops(temp_img_Kmean, 'Solidity');
% % Sort in oder of decreasing solidity.
% [sortedS, sortIndexes] = sort([measurements.Solidity], 'descend');
% % Get the solidity of the most solid blob
% highestSolidity = sortedS(1);
% % Get the label of the most solid blob
% labelWithHighestSolidity = sortIndexes(1);
% 
% size(highestSolidity)
% size(labelWithHighestSolidity)
% highestSolidity
% labelWithHighestSolidity

% bw_min = (img_Kmean == min(img_Kmean(:)));
% bw_max = (img_Kmean == max(img_Kmean(:)));
% 
% bw_min = imopen(bw_min, strel('disk', 4));
% bw_max = imopen(bw_max, strel('disk', 4));
% 
% figure(); imshow(bw_min, [])
% figure(); imshow(bw_max, [])
% 
% bw_border(bw_border == 0) = NaN;
% temp_vec_min = reshape(bw_min .* bw_border, 1, []);
% temp_vec_max = reshape(bw_max .* bw_border, 1, []);
% 
% temp_vec_min(isnan(temp_vec_min)) = [];
% temp_vec_max(isnan(temp_vec_max)) = [];


% if nnz(temp_vec_min) > nnz(temp_vec_max)
%     bw = bw_max;
% else
%     bw = bw_min;
% end

% if nnz(temp_vec == min(temp_vec(:))) > nnz(temp_vec == max(temp_vec(:)))
%     bw = (img_Kmean == max(img_Kmean(:)));
% else
%     bw = (img_Kmean == min(img_Kmean(:)));
% end

% bw_border(bw_border == 0)           = 3;
% temp_img_Kmean                      = img_Kmeans .* bw_border;
% % figure(); imshow(temp_img_Kmean, []); title('Kmeans image')
% temp_img_Kmean(temp_img_Kmean >= 3) = 3;
% temp_vec                            = reshape(temp_img_Kmean(:), 1, []);
% % temp_vec(isnan(temp_vec))           = [];
% temp_vec(temp_vec == 3)             = [];

% if nnz(temp_vec == min(temp_vec(:))) > nnz(temp_vec == max(temp_vec(:))) 
%     bw = (img_Kmeans == max(img_Kmeans(:)));
% else
%     bw = (img_Kmeans == min(img_Kmeans(:)));
% end

if nnz(bw_0) > nnz(bw_1)
    % && (nnz(img_Kmean == max(img_Kmean(:))) < 0.25*size(img,1)*size(img,2))
    bw = bw_1;
else
    bw = bw_0;
end

% figure(); imshow(bw, [])
% 
% end

%% Clear temporary variables
clear new_img bw_0 bw_1 inner_coord displacement d1 d2 temp Z Xcoord ...
    Ycoord k start_mat clusterIndexes clusterCenters img_Kmean maxValue ...
    indexOfMaxValue;

end