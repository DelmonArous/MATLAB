function [A_outer, r_outer, A_inner, r_inner] = estimateSpheroidArea( ...
    img_RGB, img_MAP, foldername, filename, sensitivity)

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
% figure(); imshowpair(RGB, img_enhanced, 'montage')
% title('Original (left) and Contrast Enhanced (right) Image')

%% Convert to grayscale image and apply threshold to create binary image
img_enhanced_grayscale = double(rgb2gray(img_enhanced));
%figure(); imshow(img__enhanced_grayscale)
img_enhanced_grayscale = adapthisteq(img_enhanced_grayscale); %, ...
%'ClipLimit', 0.02, 'Distribution', 'exponential');
%figure(); imshow(img__enhanced_grayscale)
bw_outer = imbinarize(img_enhanced_grayscale, 'adaptive', ...
    'ForegroundPolarity', 'dark', 'Sensitivity', sensitivity);

%% Plot generated binary bws
bw_outer = imclearborder(~bw_outer, 8);
figure(); imshow(bw_outer)
bw_outer = imfill(bw_outer, 'holes');
figure(); imshow(bw_outer)
bw_outer = bwmorph(bw_outer, 'bridge', Inf);
figure(); imshow(bw_outer)
bw_outer = bwmorph(bw_outer, 'fill', Inf);
figure(); imshow(bw_outer)
% % Remove objects containing fewer than 200 pixels
bw_outer = bwareaopen(bw_outer, 200);
figure(); imshow(bw_outer)
% % Remove disk-shaped structures with radius less than 85 pixels
bw_outer = imopen(bw_outer, strel('disk', 85));
figure(); imshow(bw_outer)
% Draw a perimeter around the focused structure
bw_perim_outer = bwperim(bw_outer);
overlay_outer = imoverlay(img_enhanced, bw_perim_outer, [.3 1 .3]);
h = figure(); imshow(overlay_outer)

%% Find circle of the spheroid
img_center = size(bw_outer)./2 + .5;

stats = regionprops('struct', bw_outer, 'Centroid', 'Eccentricity', ...
    'MajorAxisLength', 'MinorAxisLength');
%stats(stats.Eccentricity > 0.5, :) = [];

x_censtroid = [];
y_centroid = [];
for i = 1:length(stats)
    x_centroid(i) = stats(i).Centroid(1);
    y_centroid(i) = stats(i).Centroid(2);
end

[~, ind] = min((x_centroid - img_center(2)).^2 + ...
    (y_centroid - img_center(1)).^2);

center_outer = stats(ind).Centroid;
diameter_outer = mean([stats(ind).MajorAxisLength stats(ind).MinorAxisLength], 2);
radius_outer_px = diameter_outer ./ 2;

%% Estimate area
conversion = 1.201; % pixels/mu
r_outer = radius_outer_px .* (1/conversion);
A_outer = pi .* r_outer.^2 .* 10^(-6); % mm^2

%% Check if  necrotic center exists
if radius_outer_px > 290
    
    se = strel('disk', 110); % 250
    
    % Opening is an erosion followed by a dilation
    Io = imopen(greenChannel, se);
    %figure(); imshow(Io)
    
    % Opening-by-reconstruction is an erosion followed by a morphological 
    % reconstruction
    Ie = imerode(greenChannel, se);
    Iobr = imreconstruct(Ie, greenChannel);
    %figure(); imshow(Iobr)
    
    % Following the opening with a closing removes the dark spots and 
    Ioc = imclose(Io,se);
    %figure(); imshow(Ioc)
    
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    %figure(); imshow(Iobrcbr)
    
    % Calculate the regional minima of Iobrcbr to obtain foreground markers
    bw_inner = imregionalmin(Iobrcbr, 8);
    
    % Clean the edges of the marker blobs and then shrink them a bit by
    % closing followed by an erosion
    se2 = strel(ones(5,5));
    bw_inner = imclose(bw_inner, se2);
    bw_inner = imerode(bw_inner, se2);
    bw_inner = bwareaopen(bw_inner, 20);
    
    % Clean up some stray isolated pixels
    bw_inner = imclearborder(bw_inner, 8);
    %figure(); imshow(bw)
    bw_inner = imfill(bw_inner, 'holes');
    %figure(); imshow(bw); title('Imfill holes')
    bw_inner = bwareaopen(bw_inner, 200);
    %figure(); imshow(bw); title('Bwareaopen')
    bw_inner = bwmorph(bw_inner, 'thicken', 20);
    %figure(); imshow(bw); title('Thicken')
    bw_inner = bwmorph(bw_inner, 'bridge', Inf);
    %figure(); imshow(bw); title('Bride')
    bw_inner = imfill(bw_inner, 'holes');
    %figure(); imshow(bw); title('Imfill holes')
    bw_inner = imerode(bw_inner, strel('disk', 3));
    %figure(); imshow(bw); title('Imerode')
    bw_inner = imopen(bw_inner, strel('disk', 85));
    %figure(); imshow(bw); title('Imopen')
    
    % Estimate median intensity in the necrotic center
    ROIgreenChannelValues = reshape(greenChannel .* bw_inner, 1, []);
    ROIgreenChannelValues(ROIgreenChannelValues == 0) = [];
    ROIgreenChannelValues(isnan(ROIgreenChannelValues)) = [];
    ROIgreenChannelMedianValue = median(ROIgreenChannelValues);
    
    %[filename num2str(ROIgreenChannelMedianValue)]
    
    if ROIgreenChannelMedianValue < 0.43
        
        bw_perim_inner = bwperim(bw_inner);
        overlay_inner = imoverlay(overlay_outer, bw_perim_inner, [.3 1 .3]);
        imshow(overlay_inner)
        
        stats = regionprops('struct', bw_inner, 'Centroid', 'Eccentricity', ...
            'MajorAxisLength', 'MinorAxisLength');
        
        center_inner = stats.Centroid;
        diameter_inner = mean([stats.MajorAxisLength stats.MinorAxisLength], 2);
        radius_inner_px = diameter_inner ./ 2;
        
        hold on
        viscircles(center_outer, radius_outer_px, 'EdgeColor', 'b', 'LineWidth', 0.5);
        plot(center_outer(1), center_outer(2), 'b*')
        plot(linspace(center_outer(1), center_outer(1)+radius_outer_px, 1000), ...
            linspace(center_outer(2), center_outer(2), 1000), 'b', ...
            'LineWidth', 0.5)
        
        viscircles(center_inner, radius_inner_px, 'EdgeColor', 'r', ...
            'LineWidth', 0.5);
        plot(center_inner(1), center_inner(2), 'r*')
        plot(linspace(center_inner(1), center_inner(1), 1000), ...
            linspace(center_inner(2), center_inner(2)+radius_inner_px, 1000), ...
            'r', 'LineWidth', 0.5)
        hold off
        
        r_inner = radius_inner_px .* (1/conversion);
        A_inner = pi .* r_inner.^2 .* 10^(-6); % mm^2
        
        str = {'Spheroid dimensions estimate', ...
            ['Growth: ' foldername(1:2) ' days'], ...
            ['Filename: ' filename], ...
            ['A_{outer} = ' num2str(round(A_outer, 4)) ...
            ' mm^2, r_{outer} = ' num2str(round(r_outer, 2)) ' \mum'], ...
            ['A_{inner} = ' num2str(round(A_inner, 4)) ...
            ' mm^2, r_{inner} = ' num2str(round(r_inner, 2)) ' \mum']};
        str{1} = ['\bf ', str{1}, ' \rm'];
        a = annotation('textbox', [0.1 0.4 0.55 0.55], 'String', str, 'FitBoxToText', 'on');
        a.FontSize = 11;
        
    else
        
        hold on
        viscircles(center_outer, radius_outer_px, 'EdgeColor', 'b', 'LineWidth', 0.5);
        plot(center_outer(1), center_outer(2), 'b*')
        plot(linspace(center_outer(1),center_outer(1)+radius_outer_px,1000), ...
            linspace(center_outer(2),center_outer(2),1000), 'b', 'LineWidth', 0.5)
        hold off
        
        r_inner = 0;
        A_inner = 0;
        
        str = {'Spheroid dimensions estimate', ...
            ['Growth: ' foldername(1:2) ' days'], ...
            ['Filename: ' filename], ...
            ['A_{outer} = ' num2str(round(A_outer, 4)) ...
            ' mm^2, r_{outer} = ' num2str(round(r_outer, 2)) ' \mum']};
        str{1} = ['\bf ', str{1}, ' \rm'];
        a = annotation('textbox', [0.1 0.4 0.55 0.55], 'String', str, 'FitBoxToText', 'on');
        a.FontSize = 11;
                
    end
    
else
    
    hold on
    viscircles(center_outer, radius_outer_px, 'EdgeColor', 'b', 'LineWidth', 0.5);
    plot(center_outer(1), center_outer(2), 'b*')
    plot(linspace(center_outer(1),center_outer(1)+radius_outer_px,1000), ...
        linspace(center_outer(2),center_outer(2),1000), 'b', 'LineWidth', 0.5)
    hold off
    
    r_inner = 0;
    A_inner = 0;
    
    str = {'Spheroid dimensions estimate', ...
        ['Growth: ' foldername(1:2) ' days'], ...
        ['Filename: ' filename], ...
        ['A_{outer} = ' num2str(round(A_outer, 4)) ...
        ' mm^2, r_{outer} = ' num2str(round(r_outer, 2)) ' \mum']};
    str{1} = ['\bf ', str{1}, ' \rm'];
    a = annotation('textbox', [0.1 0.4 0.55 0.55], 'String', str, 'FitBoxToText', 'on');
    a.FontSize = 11;
    
end

% saveas(h, sprintf('%s_%s.jpg', foldername, filename))

end