function [bw] = borderExtractor(img_norm)

figure(); imshow(img_norm, [])
hsv = rgb2hsv(img_norm);
v = hsv(:,:,3);

%% Preprocessing of the image
A = adapthisteq(img_norm, 'NumTiles', [16 16], 'ClipLimit', 0.02, ...
    'Distribution', 'exponential');
% figure(); imshow(A, [])

if ~isempty(size(img_norm, 3)) && size(img_norm,3) == 3
    im = double(mean(A,3));
elseif ~isempty(size(img_norm,3)) && size(img_norm,3) == 2
    im = double(mean(A,2));
else
    im = double(A);
end

% figure(); imshow(im, [])
% title('im')
filled = imfill(im, 'holes');
% figure(); imshow(filled, [])
% title('filled')

%% Threshold detection using most occuring value in histogram
[counts, x] = hist(img_norm(:), 255);
bw1 = img_norm <= x(find(counts == max(max(counts)))) + 0.2;
% figure(); imshow(bw1, [])
bw2 = img_norm >= x(find(counts == max(max(counts)))) - 0.2;
% figure(); imshow(bw2, [])
bw = bw1 & bw2;
bw = imclearborder(bw, 8);
%figure(); imshow(bw, [])
% title('bw')

% bw(1,:) = 0; %bw(:,1) = 0;
% bw(end,:) = 0; %bw(:,end) = 0;

% Clean up
bw = bwmorph(bw, 'bridge', Inf);
bw = imclearborder(bw, 8);
bw = bwmorph(bw, 'clean', Inf);
bw = imopen(bw, strel('disk', 10));
bw = bwareaopen(bw, 30000);
bw = imclose(bw, strel('disk', 200));
bw = imfill(bw, 'holes');
% figure(); imshow(bw, [])
% title('bw')

%% 
% bw_border = imbinarize(img, graythresh(img)); % ;
% figure(); imshow(bw_border, [])
% 
% % Remove pixels along border
% out = imclearborder(bw_border);
% figure(); imshow(out, []);
% 
% % Obtain pixels that are along border
% bw_border2 = bw_border; % make copy
% bw_border2(out) = 0; % set pixels not belonging to boundary to 0
% figure(); imshow(bw_border2, []);
% 
% % Fill holes for regions
% % out = bwmorph(out, 'bridge', Inf);
% % out_fill = imfill(out, 'holes');
% % out_fill = bwmorph(out_fill, 'clean', Inf);
% out2_fill = ~bwareaopen(~bw_border2, 500);
% figure(); imshow(out2_fill, [])
% out2_fill = bwmorph(out2_fill, 'majority', Inf);
% figure(); imshow(out2_fill, [])
% out2_fill = imfill(out2_fill, 'holes');
% figure(); imshow(out2_fill, [])
% bw_border = ~out2_fill;
% 
% % out2_fill = bwareafilt(out2_fill, 1);
% % figure(); imshow(out2_fill, []);
% % title('Test4')
% 
% % bw_border(1,:) = 0; bw_border(:,1) = 0;
% bw_border(end,:) = 0; % bw_border(:,end) = 0;
% bw_border = imclearborder(bw_border, 8);
% bw_border = imclose(bw_border, strel('disk', 20));
% bw_border = imfill(bw_border, 'holes');
% bw = bw_border;
% figure(); imshow(bw, [])

%% Clean up
% bw = bwmorph(bw, 'bridge', Inf);
% bw = imclearborder(bw, 8);
% bw = bwmorph(bw, 'clean', Inf);
% bw = imfill(bw, 'holes');
% bw = bwmorph(bw, 'majority', Inf);
% bw = imopen(bw, strel('disk', 5));
% bw = bwareaopen(bw, 1500);
% figure(); imshow(bw, [])
% title('bw')

% bw = imcomplement(bw);
% bw = bwmorph(bw, 'bridge', Inf);
% bw = imopen(bw, strel('disk', 3));
% bw = bwareaopen(bw, 450);
% bw = imcomplement(bw);
% figure(); imshow(bw, [])

end