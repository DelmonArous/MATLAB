function computeCuvetteIntensityMap(filename_img, filename_intensity, ...
    px_crop_start, px_crop_end, px_surf_start, px_surf_end)

conversion = 1/49.0196; % in mm/pixel

% Filenames
[~, img_name, ~] = fileparts(filename_img);
[~, colormap_name, ~] = fileparts(filename_intensity);

% Read image, get dimensionf of the image, and extract the individual
% red, green, and blue color channels.
img = flip(imread(filename_img), 2);
img = imcrop(img, [px_crop_start(1) px_crop_start(2) ...
    px_crop_end(1)-px_crop_start(1) px_crop_end(2)-px_crop_start(2)]);
[rows, columns, numberOfColorChannels] = size(img);
redChannel = img(:,:,1);
greenChannel = img(:,:,2);
blueChannel = img(:,:,3);

% Read the intensity image
SI_img = flip(imread(filename_intensity),2);
SI_img = imcrop(SI_img, [px_crop_start(1) px_crop_start(2) ...
    px_crop_end(1)-px_crop_start(1) px_crop_end(2)-px_crop_start(2)]);
SI_img = double(rgb2gray(SI_img));
% SI_img = double(255*mat2gray(rgb2gray(SI_img)));
% imadjust(rgb2gray(SI_img), [0 18], [0 1])
% SI_img = double(imadjust(rgb2gray(SI_img), [0 18], [0 255]));
% SI_img(isnan(SI_img)) = 0; SI_img(isinf(SI_img)) = 0;
% SI_img = log(SI_img);

% SI_img(1019:1049, 1:9) = 0;
% SI_img(1190:1213, 1:3) = 0;

% Create a spatial referencing object associated with the image, and
% use the referencing object to set the x- and y-axes limits in the
% world coordinate system
sizex = size(img, 2);
sizey = size(img, 1);
xmax = sizex * conversion;
ymax = sizey * conversion;
RI = imref2d(size(img));
RI.XWorldLimits = [0 xmax];
RI.YWorldLimits = [0 ymax];

% Get mask of the intensity image and project it onto the color image
mask = SI_img > 0;
SI_img_mask = SI_img .* mask;
%SI_img_mask = smooth2a(SI_img_mask, 1, 1);
maskedRgbImage = bsxfun(@times, img, cast(~mask, 'like', img));

%% Surface plot
h = figure();
s = surf(SI_img .* (SI_img > 0));
xData = s.XData;
yData = s.YData;
zData = s.ZData;
xData = xData(xData(1):xData(125));
yData = yData(px_surf_start:px_surf_end);
zData = zData(px_surf_start:px_surf_end, xData(1):xData(125));
zData(isnan(zData)) = 0;

% for i = 1:size(zData,1) % 1100:1113 %
%     
%     SI = zData(i,:);
%     diffSI = abs(diff(SI));
%     thres = quantile(diffSI, 0.25);
%     tempmax = max(SI(:));
%     tempmin = min(SI(:));
%     
% %     figure();
% %     plot(xData, SI, '-')
% %     xlabel('Wavelength')
% %     ylabel('Signal Intensity')
% %     ylim([tempmin tempmax])
%     
% %     for k = 1
% %         
% %         for j = (1 + k):(length(SI) - k)
% %             if abs(SI(j) - SI(j+k)) > thres
% %                 SI(j) = NaN;
% %                 SI(j+k) = NaN;
% %             elseif abs(SI(j) - SI(j-k)) > thres
% %                 SI(j) = NaN;
% %                 SI(j-k) = NaN;
% %             end
% %         end
% %         SI = fillmissing(SI, 'pchip');
% %                 
% %     end
%     
%     SI = fillmissing(SI, 'pchip');
%     windowWidth = 3;
%     SI = movmin(SI, windowWidth); % movemean
%     windowWidthSG = 7;
%     SI = sgolayfilt(SI, 3, windowWidthSG);
% %      figure();
% %     plot(xData, SI, '-')
% %     xlabel('Wavelength')
% %     ylabel('Signal Intensity')    
% %     ylim([tempmin tempmax])
%     
%     %SI = imopen(SI, ones(25,1));
%     zData(i,:) = SI;
%     
%     %zData(i,:) = smooth(xData, zData(i,:), 'sgolay', 1);
%     % new_zData(i,:) = pchip(xData, zData(i,:), new_xData);
%     %new_zData(i,:) = smooth(new_xData, new_zData(i,:));
%     %new_zData(i,:) = medfilt1(new_zData(i,:), 10, 'omitnan','truncate');
% end

% % new_xData = xData(1):0.1:xData(end); % linspace(xData(1), xData(end), 10*length(xData)); %
% % new_yData = (yData(1):0.1:yData(end)).'; % linspace(yData(1), yData(end), 10*length(yData)).'; %
% % new_zData = griddata(xData, yData, zData, new_xData, new_yData, 'cubic');
% % %new_zData = interp2(xData, yData, zData, new_xData, new_yData, 'spline');
% % xData = new_xData; yData = new_yData; zData = new_zData;
%
% % %zData = medfilt2(zData, [5 5]);
% % % zData = smoothn(zData);
% % % zData = smooth2a(zData, 10, 10);
% %
% % % k = 1;
% % % verticallySmoothedImage = sgolayfilt(zData, 1, 21, [], 1);
% % % zData = sgolayfilt(verticallySmoothedImage, 1, 21, [], 2);
% %
% % % g = gausswin(20);
% % % g = g/sum(g);
% % % for i = 1:size(zData,1)
% % %     %zData(i,:) = conv(zData(i,:), g, 'same');
% % %     %zData(i,:) = smooth(xData, zData(i,:), 'sgolay', 1);
% % %     new_zData(i,:) = pchip(xData, zData(i,:), new_xData);
% % %     %new_zData(i,:) = smooth(new_xData, new_zData(i,:));
% % %     %new_zData(i,:) = medfilt1(new_zData(i,:), 10, 'omitnan','truncate');
% % % end
% % % for i = 1:size(zData,2)
% % %     %zData(:,i) = conv(zData(:,i), g, 'same');
% % %     %zData(:,i) = smooth(yData, zData(:,i), 'sgolay', 1);
% % %     new_zData(:,i) = pchip(xData, zData(:,i), new_xData);
% % % end
% %
% % % new_zData = gridfit(xData, yData.', zData, new_xData, new_yData);
% % % new_zData = interp2(xData, yData, zData, new_xData, new_yData, 'cubic');
% % % xData = new_xData; yData = new_yData; zData = new_zData;

h = figure();
set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
s = surf(xData, yData, zData, 'LineStyle', 'none', 'FaceColor', 'interp');
xlim([xData(1) xData(end)])
ylim([yData(1) yData(end)])
zlim([0 ceil(max(zData(:)))]); % ceil(max(zData(:)))
colormap(jet(4096))
xlabel('mm'); ylabel('mm'); zlabel('Intensity');
set(gca,'Ydir', 'reverse')
s.EdgeColor = 'none';
grid on
shading interp
view(3); axis vis3d; % camlight HEADLIGHT;

% Change x and y tick labels from pixels to mm
% xticks = get(gca, 'XTick');
% set(gca,'XTickLabel', num2str(round(conversion.*xticks',1)))
addMM_x = @(x) sprintf('%.1f', x * conversion);
addMM_y = @(y) sprintf('%.0f', y * conversion);
xticklabels(cellfun(addMM_x,num2cell(xticks'), 'UniformOutput', false));
yticklabels(cellfun(addMM_y,num2cell(yticks'), 'UniformOutput', false));

% Find maxium intensity (Bragg peak) region and plot it
hold on
maxValue = max(zData(:));
[rowsOfMaxes, colsOfMaxes] = find(zData == maxValue);
Bragg_mask = zData .* (zData == maxValue);
y_minind = min(rowsOfMaxes(:));
y_maxind = max(rowsOfMaxes(:));
y_medianind = round(median(rowsOfMaxes(:)));
x_minind = min(colsOfMaxes(:));
x_maxind = max(colsOfMaxes(:));
x_medianind = round(median(colsOfMaxes(:)));
[x_minind x_maxind round(median(colsOfMaxes(:))) max(zData(:)) ...
    round(median(colsOfMaxes(:)))*conversion]

% % p = zeros(1,3);
% % p(1) = plot3(linspace(xData(x_minind), xData(x_minind), ...
% %     length(yData)), yData, zData(:, x_minind), 'Linestyle', '-.', ...
% %     'Color', 'black', 'LineWidth', 2.0, 'DisplayName', 'Estimated Bragg peak region');
% % p(2) = plot3(linspace(xData(x_maxind), xData(x_maxind), ...
% %     length(yData)), yData, zData(:, x_maxind), 'Linestyle', '-.', ...
% %     'Color', 'black', 'LineWidth', 2.0, 'DisplayName', 'Estimated Bragg peak region');

plot3(linspace(xData(x_medianind), xData(x_medianind), ...
    length(yData)), yData, zData(:, x_medianind), 'Linestyle', '-.', ...
    'Color', 'black', 'LineWidth', 2.0, 'DisplayName', 'Estimated Bragg peak region');
% contour3(xData, yData, zData, [maxValue maxValue], 'Linestyle', '-.', ...
%     'Color', 'black', 'LineWidth', 2.0, 'DisplayName', 'Estimated Bragg peak region');
% legend(p(1), 'Estimated width of maximum intensity', 'Location', 'North')
hold off
set(gca, 'FontSize', 20)
% saveas(h, sprintf('SurfImageTop_%s_on_%s.jpg', ...
%     colormap_name(1:2), img_name(1:2)));

zData(y_medianind, :)

%% Impose colormap on image and plot
h = figure();
set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
imshow(maskedRgbImage, RI)
hold on
ih = imshow(SI_img_mask, RI, [], 'colormap', jet(4096));
set(ih, 'AlphaData', SI_img_mask > 0);
colormap(jet(4096))
c = colorbar;
c.Label.String = 'Intensity';
caxis([0 ceil(max(zData(:)))]) %
xlabel('mm')
ylabel('mm')
shading interp
set(gca, 'FontSize', 14)
% saveas(h, sprintf('ColormapImage_%s_on_%s.jpg', ...
%     colormap_name(1:2), img_name(1:2)));

end