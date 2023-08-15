function bw = extractBLOBs(img, border, param)

% The following variables are required for execution: 
%   img = input image that is either a grayscale or a single color (red) 
%       channel image of the inquired cell colonies. 
%   border = binary mask of the cell container boundaries.
%   area_min = minimum user-defined colony area. 
%
% The following variables are returned upon successful completion when 
% input arguments are provided:
%   bw = binary mask containing the BLOBs representing the cell colony 
%       conglomerations. 

%% Min-max Normalization between 0-1

% % temp_img = img;
% % img = img .* border;
% % img = (img - min(min(img))) ./ abs(max(max(img)) - min(min(img)));
% % img = (adapthisteq(img) - min(min(adapthisteq(img)))) ./ ...
% %     abs(max(max(adapthisteq(img))) - min(min(adapthisteq(img))));
% % figure(); imshow(img, [])

% img = (img - min(min(img))) ./ abs(max(max(img)) - min(min(img)));

% img1 = img .* border;

% % Opening-by-Reconstruction
% img_e       = imerode(img, strel('disk', param.obrcbr_size_pca));
% img_obr     = imreconstruct(img_e, img);
% 
% % Opening-Closing-by-Reconstruction
% img_obrd    = imdilate(img_obr, strel('disk', param.obrcbr_size_pca));
% img_obrcbr  = imreconstruct(imcomplement(img_obrd), imcomplement(img_obr));
% bckg        = imcomplement(img_obrcbr);
% 
% img = img - bckg;
% img = imgaussfilt(img, param.gaussfilt_size_pca); 
% img = img .* border; 
% % img = (adapthisteq(img) - min(min(adapthisteq(img)))) ./ ...
% %     abs(max(max(adapthisteq(img))) - min(min(adapthisteq(img))));

% bckg1 = medfilt2(img, [70 70], 'symmetric');
% figure(); imshow(bckg1, [])
% % bckg2 = ordfilt2(img, round(100^2/2), true(100));
% % figure(); imshow(bckg2, [])
% % figure(); imshow(img, [])
% % img = img - bckg2;
% % figure(); imshow(img, [])

% figure(); imshow(img, [])
img = adapthisteq(img);
img = ( img - min(img(:)) ) ./ abs( max(img(:)) - min(img(:)) );

%% Threshold by means of K-means to extract conglomerate colony BLOBs 
bw = runKmeans(img);
bw = bw .* border; % added!!!!

%% Clean up by applying morphological operations

if param.intensiveseg == 1
    a_max = 100000000000000;
else
    a_max = 2*param.area(4);
end

bw = imclearborder(bw, 8);
bw = bwmorph(bw, 'clean', Inf);
bw = bwmorph(bw, 'hbreak', Inf);
bw = bwmorph(bw, 'spur', Inf);
bw = imclose(bw, strel('disk', 3));
bw = imerode(bw, strel('disk', 2));
bw = bwpropfilt(bw, 'Area', [param.area(1)/2 a_max]);
% bw = bwmorph(bw, 'bridge', Inf);
% bw = imfill(bw, 'holes');
% bw = bwpropfilt(bw, 'Area', [param.area(1) 10000000]);
% bw = imfill(bw, 'holes');
bw = imdilate(bw, strel('disk', 3));
bw = bwmorph(bw, 'majority');
bw = bwmorph(bw, 'thicken');

% Remove objects containing fewer than a given pixels
% bw = bwareaopen(bw, param.area(1)); % increase to remove more "smudges"
bw = bwpropfilt(bw, 'Area', [param.area(1)/2 a_max]);
% bw = bwpropfilt(bw, 'Eccentricity', [0 0.99]);

% if strcmp(param.pca, 'PCA2')
%     bw = conv2(single(bw), ones(11) / 11^2, 'same');
%     bw = bw > 0.5;
% end

end