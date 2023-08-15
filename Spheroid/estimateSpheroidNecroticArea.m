function [area, center, radii_um] = estimateSpheroidNecroticArea( ...
    img_RGB, img_MAP, foldername, filename)

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
greenChannel = img_enhanced(:, :, 2);

%%
se = strel('disk', 250);
Io = imopen(greenChannel, se);
%figure(); imshow(Io)

Ie = imerode(greenChannel, se);
Iobr = imreconstruct(Ie, greenChannel);
%figure(); imshow(Iobr)

Ioc = imclose(Io,se);
%figure(); imshow(Ioc)
 
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
%figure(); imshow(Iobrcbr)
 
mask = imregionalmin(Iobrcbr, 8);

I2 = labeloverlay(greenChannel, mask);
%figure(); imshow(I2)
 
se2 = strel(ones(5,5));
mask = imclose(mask, se2);
mask = imerode(mask, se2);
mask = bwareaopen(mask, 20);
I3 = labeloverlay(greenChannel, mask);
%figure(); imshow(I3)

mask = imclearborder(mask, 8);
%figure(); imshow(mask)
mask = imfill(mask, 'holes');
%figure(); imshow(mask); title('Imfill holes')
mask = bwareaopen(mask, 200);
%figure(); imshow(mask); title('Bwareaopen')
mask = bwmorph(mask, 'thicken', 15);
%figure(); imshow(mask); title('Thicken')
mask = bwmorph(mask, 'bridge', Inf);
%figure(); imshow(mask); title('Bride')
mask = imfill(mask, 'holes');
%figure(); imshow(mask); title('Imfill holes')
mask = imerode(mask, strel('disk', 3));
%figure(); imshow(mask); title('Imerode')
mask = imopen(mask, strel('disk', 85)); 
%figure(); imshow(mask); title('Imopen')

mask_perim = bwperim(mask);
overlay = imoverlay(img_enhanced, mask_perim, [.3 1 .3]);
h = figure(); imshow(overlay)

%%
stats = regionprops('struct', mask, 'Centroid', 'Eccentricity', ...
    'MajorAxisLength', 'MinorAxisLength');

center = stats.Centroid;
diameter = mean([stats.MajorAxisLength stats.MinorAxisLength], 2);
radii = diameter ./ 2;

hold on
viscircles(center, radii, 'EdgeColor', 'b', 'LineWidth', 0.5);
plot(center(1), center(2), 'b*')
plot(linspace(center(1),center(1)+radii,1000), ...
    linspace(center(2),center(2),1000), 'b', 'LineWidth', 0.5)
hold off

%% Estimate area
conversion = 1.201; % pixels/mu 
radii_um = radii .* (1/conversion);
area = pi .* radii_um.^2 .* 10^(-6); % mm^2 

%%
str = {'Spheroid dimensions estimate', ...
    ['Growth: ' foldername(1:2) ' days'], ...
    ['Filename: ' filename], ...
    ['A = ' num2str(round(area, 4)) ' mm^2, r = ' num2str(round(radii_um, 2)) ' \mum']};
str{1} = ['\bf ', str{1}, ' \rm'];
a = annotation('textbox', [0.1 0.4 0.5 0.5], 'String', str, 'FitBoxToText', 'on');
a.FontSize = 12;

%% Compute the gradient magnitude as the Segmentation Function
% gmag = imgradient(img__enhanced_grayscale);
% figure(2); imshow(gmag, [])

% The image can not be segmented by using the watershed transform directly 
% on the gradient magnitude. This results in oversegmentation
% L = watershed(gmag);
% Lrgb = label2rgb(L);
% figure(); imshow(Lrgb) 

%% Mark the forground objects
% se = strel('disk', 20); % 20

% % Remove disk-shaped structures with radius less than 85 pixels
% Io = imopen(img_enhanced_grayscale, se);
% figure(3); imshow(Io)
%  
% % Compute the opening-by-reconstruction using imerode and imreconstruct
% Ie = imerode(img_enhanced_grayscale, se);
% Iobr = imreconstruct(Ie, img_enhanced_grayscale);
% figure(4); imshow(Iobr)
% 
% % Closing to remove dark spots and stem marks
% Ioc = imclose(Io, se);
% figure(5); imshow(Ioc)
% 
% Iobrd = imdilate(Iobr, se);
% Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
% Iobrcbr = imcomplement(Iobrcbr);
% figure(6); imshow(Iobrcbr)
% 
% % Calculate the regional maxima of Iobrcbr to obtain foreground markers
% fgm = imregionalmax(Iobrcbr);
% figure(7); imshow(fgm)
% 
% % Superimpose the foreground marker image on the original image
% I2 = labeloverlay(img_enhanced_grayscale, fgm);
% figure(8); imshow(I2)
% 
% % Clean the edges of the marker blobs and then shrink them a bit
% se2 = strel(ones(5,5));
% fgm2 = imclose(fgm, se2);
% fgm3 = imerode(fgm2, se2);
% 
% % Remove all blobs that have fewer than a certain number of pixels 
% fgm4 = bwareaopen(fgm3, 20);
% figure(9); imshow(fgm4);
% 
% mask_perim = bwperim(fgm4);
% overlay = imoverlay(img_enhanced, mask_perim, [.3 1 .3]);
% h = figure(); imshow(overlay)
% 
% [centers, radii] = imfindcircles(fgm4, [200 350], ...
%     'ObjectPolarity', 'dark', 'Sensitivity', 0.92, 'Method', 'TwoStage'); % 
% 
% hold on
% viscircles(centers, radii, 'Color', 'b');
% hold off

% stats = regionprops('table', fgm4, 'Centroid', 'Eccentricity', ...
%     'MajorAxisLength', 'MinorAxisLength');
% 
% center = stats.Centroid;
% diameter = mean([stats.MajorAxisLength stats.MinorAxisLength], 2);
% radii = diameter ./ 2;
% 
% hold on
% viscircles(center, radii, 'EdgeColor', 'b', 'LineWidth', 0.5);
% hold off

% %% Compute background markers
% mask = imbinarize(Iobrcbr);
% figure(11); imshow(~mask)
% 
% D = bwdist(~mask);
% DL = watershed(D);
% bgm = DL == 0;
% figure(12); imshow(bgm)
% 
% %% Compute the Watershed Transform of the Segmentation Function
% % Modify the gradient magnitude image so that its only regional minima
% % occur at foreground and background marker pixels.
% gmag2 = imimposemin(gmag, bgm | fgm4);
% % Compute the watershed-based segmentation
% L = watershed(gmag2);
% 
% %% Visualize the Result
% labels = imdilate(L == 0, ones(3,3)) + 2*bgm + 3*fgm4;
% I4 = labeloverlay(img__enhanced_grayscale, labels);
% figure(13); imshow(I4)
% 
% Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% figure(); imshow(Lrgb)
% 
% figure();
% imshow(img__enhanced_grayscale)
% hold on
% himage = imshow(Lrgb);
% himage.AlphaData = 0.3;

%% Return the top-hat and bottom-hat transforms, respectively
% Itop = imtophat(img__enhanced_grayscale, se);
% Ibot = imbothat(img__enhanced_grayscale, se);
% 
% Ienhance = imsubtract(imadd(Itop, img__enhanced_grayscale), Ibot);
% figure(); imshow(Ienhance)
% 
% Iec = imcomplement(Ienhance);
% figure(); imshow(Iec)
% 
% Iemin = imextendedmin(Iec, 1);
% Iimpose = imimposemin(Iec, ~Iemin);
% figure(); imshow(~Iemin)
% 
% wat = watershed(Iimpose);
% rgb2 = label2rgb(wat);
% figure(); imshow(rgb2);

%%
% mask = imbinarize(img_enhanced_grayscale, 'adaptive', ...
%      'ForegroundPolarity', 'dark', 'Sensitivity', 0.57); 
% % %mask = im2bw(img_enhanced_grayscale, graythresh(img_enhanced_grayscale));
% figure(); imshow(mask)
%  
% mask = imclearborder(mask, 8);
% figure(); imshow(mask)
% mask = bwmorph(mask, 'bridge', Inf);
% figure(); imshow(mask); title('1st Morph Bridge')
% mask = imfill(mask, 8, 'holes');
% figure(); imshow(mask); title('Imfill')
% mask = imopen(mask, strel('disk', 3)); 
% figure(); imshow(mask); title('Imopen')
% mask = imclose(mask, strel('disk', 5));
% figure(); imshow(mask); title('Imclose')
% mask = imfill(mask, 8, 'holes');
% figure(); imshow(mask); title('Imfill') 
% % mask = bwmorph(mask, 'majority', Inf);
% % figure(); imshow(mask); title('Majority')
% mask = bwmorph(mask, 'thicken', 3);
% figure(); imshow(mask); title('Thicken')
% mask = bwmorph(mask, 'bridge', Inf);
% figure(); imshow(mask); title('2nd Morph Bridge')
% mask_perim = bwperim(mask);
% overlay1 = imoverlay(img_enhanced_grayscale, mask_perim, [.3 1 .3]);
% figure(); imshow(overlay1)

% mask_em = imextendedmin(img_enhanced_grayscale, 0.3);
% figure(); imshow(mask_em)
% mask_em = imclose(~mask_em, strel('disk', 5));
% figure(); imshow(mask_em); title('Imclose')
% mask_em = imerode(mask_em, strel('disk', 5));
% figure(); imshow(mask_em); title('Imerode')
% mask_em = imopen(mask_em, ones(5,5)); 
% figure(); imshow(mask_em); title('Imopen')
% mask_em = bwmorph(mask_em, 'hbreak', Inf);
% figure(); imshow(mask_em); title('Bwmorph hbreak')
% mask_em = bwmorph(mask_em, 'fill', Inf);
% figure(); imshow(mask_em); title('Bwmorph fill')
% mask_em = bwareaopen(mask_em, 85);
% figure(); imshow(mask_em); title('Bwareaopen')
% overlay2 = imoverlay(img_enhanced_grayscale, mask_perim | mask_em, [.3 1 .3]);
% figure(); imshow(overlay2)
 
% I_eq_c = imcomplement(img_enhanced_grayscale);
% I_mod = imimposemin(I_eq_c, ~mask | mask_em);
% L = watershed(I_mod);
% figure(); imshow(label2rgb(L))

% stats = regionprops('table', mask, 'Centroid', 'Eccentricity', ...
%     'MajorAxisLength', 'MinorAxisLength');
% center = stats.Centroid;
% diameter = mean([stats.MajorAxisLength stats.MinorAxisLength], 2);
% radii = diameter ./ 2;
% 
% hold on
% viscircles(center, radii, 'EdgeColor', 'b', 'LineWidth', 0.5);
% hold off 

end