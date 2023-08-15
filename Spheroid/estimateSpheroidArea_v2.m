function A_outer = estimateSpheroidArea_v2( ...
    img_RGB, img_MAP, foldername, filename, value)

%% Read image
[X, map] = rgb2ind(img_RGB, img_MAP);

% Convert the indexed image into a truecolor (RGB) image,
% then convert the RGB image into the L*a*b* color space
RGB = ind2rgb(X, map);
LAB = rgb2lab(RGB);

% Scale values to the range expected by the adapthisteq function, [0 1]
max_luminosity = 100;
L = LAB(:,:,1)/max_luminosity;

% Perform CLAHE on the L channel. Scale the result to get back
% to the range used by the L*a*b* color space
LAB(:,:,1) = adapthisteq(L, 'NumTiles', [64 64], 'ClipLimit', 0.005) ...
    * max_luminosity;
img_enhanced = lab2rgb(LAB);

% Extract the green channel from the enhanced image
greenChannel = img_enhanced(:,:,2);
%figure(); imshow(greenChannel)

%% Plot original and enhanced image
% figure(); imshowpair(RGB, LAB, 'montage')
% title('Original (left) and LAB (right) Image')
%figure(); imshowpair(RGB, img_enhanced, 'montage')
% title('Original (left) and Contrast Enhanced (right) Image')

%% Convert to grayscale image and apply threshold to create binary image
img_enhanced_grayscale = double(rgb2gray(img_enhanced));
%figure(); imshow(img__enhanced_grayscale)
img_enhanced_grayscale = adapthisteq(img_enhanced_grayscale); %, ...
%'ClipLimit', 0.02, 'Distribution', 'exponential');
%figure(); imshow(img__enhanced_grayscale)
bw_outer = zeros(size(img_enhanced_grayscale,1), size(img_enhanced_grayscale,2));
sensitivity = 0.68;

while all(bw_outer(:) == 0) || nnz(bw_outer(:)) < 20000 
        
    %sensitivity
    
    bw_outer = imbinarize(img_enhanced_grayscale, 'adaptive', ...
        'ForegroundPolarity', 'dark', 'Sensitivity', sensitivity);
    
    %% Plot generated binary bw's
    bw_outer = imclearborder(~bw_outer, 8);
    %figure(); imshow(bw_outer)
    bw_outer = bwmorph(bw_outer, 'clean', Inf);
    %figure(); imshow(bw_outer)
    bw_outer = imopen(bw_outer, strel('disk', 2));
    %figure(); imshow(bw_outer)
    bw_outer = imfill(bw_outer, 'holes');
    %figure(); imshow(bw_outer)
    bw_outer = bwmorph(bw_outer, 'bridge', Inf);
    %figure(); imshow(bw_outer)
    bw_outer = bwmorph(bw_outer, 'fill', Inf);
    %figure(); imshow(bw_outer)
    bw_outer = imclose(bw_outer, strel('disk', 10));
    %figure(); imshow(bw_outer)
    bw_outer = imfill(bw_outer, 'holes');    
    %figure(); imshow(bw_outer)
    % % Remove objects containing fewer than 15000 pixels
    bw_outer = bwareaopen(bw_outer, 15000);
    %figure(); imshow(bw_outer)
    % % Remove disk-shaped structures with radius less than 'value' pixels
    bw_outer = imopen(bw_outer, strel('disk', 75));
    %figure(); imshow(bw_outer)
    
    sensitivity = sensitivity - 0.01;

end

%nnz(bw_outer(:))
%figure(); imshow(bw_outer)

%% Draw a perimeter around the final focused structure
bw_perim_outer = bwperim(bw_outer);
overlay_outer = imoverlay(img_enhanced, bw_perim_outer, [.3 1 .3]);
h = figure(); imshow(overlay_outer);

%% Estimate area
conversion = 1.201; % pixels/mu
%r_outer = radius_outer_px .* (1/conversion);
%A_outer = pi .* r_outer.^2 .* 10^(-6); % mm^2
A_outer = nnz(bw_outer) ./ conversion^2 .* 10^(-6); % mm^2 

str = {'Spheroid dimensions estimate', ...
    ['Growth: ' foldername(8:9) ' days'], ...
    ['Filename: ' filename], ...
    ['A_{outer} = ' num2str(round(A_outer, 4)) ' mm^2']};
str{1} = ['\bf ', str{1}, ' \rm'];
a = annotation('textbox', [0.1 0.4 0.55 0.55], 'String', str, 'FitBoxToText', 'on');
a.FontSize = 11;
    
saveas(h, sprintf('%s_%s.jpg', foldername, filename))

end