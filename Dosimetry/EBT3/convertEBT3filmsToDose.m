function [struct] = convertEBT3filmsToDose(path, calib, channel, ROI_bb) 

%% Find and store all image files in directory
filelist = getAllFiles(path);

%%

for i = 1:length(filelist)
    
    [~, fn, e] = fileparts(filelist{i});
        
    struct{i}.filename = [fn e];
    struct{i}.(channel) = channel;
    
    % Read EBT3 film images
    [~, I_red, I_green, I_blue, I_gray] = readEBT3film(filelist{i});
        
    if strcmp(channel, 'red')
        struct{i}.img = I_red;
    elseif strcmp(channel, 'green')
        struct{i}.img = I_green;
    elseif strcmp(channel, 'blue')
        struct{i}.img = I_blue;
    else
        struct{i}.img = I_gray;
    end
    
    % Convert pixel values into netOD, and then netOD into dose
    struct{i}.img_netOD = calculateNetOD(calib.PV_control, ...
        struct{i}.img, calib.PV_bckg);
    struct{i}.img_dose = model1(struct{i}.img_netOD, ...
        calib.mdl.Coefficients.Estimate(1), ...
        calib.mdl.Coefficients.Estimate(2), calib.n);
    
    struct{i}.ROI = imcrop(struct{i}.img_dose, [ROI_bb{i} 100 100]);
    
%     bw = imbinarize(struct{i}.img_netOD);
%     if i == 7 || i == 8 || i == 9 || i == 10
%         bw = bwmorph(bw, 'clean', Inf);
%         bw = bwmorph(bw, 'hbreak', Inf);
%         bw = bwmorph(bw, 'spur', Inf);
%         bw = imclose(bw, strel('disk', 4));
%         bw = imerode(bw, strel('disk', 2));
%         bw = bwmorph(bw, 'bridge', Inf);
%         bw = imdilate(bw, strel('disk', 3));
%         bw = bwpropfilt(bw, 'Eccentricity', [0 0.975]);
%         bw = imopen(bw, strel('disk', 20));
%         bw = imfill(bw, 'holes');
%     else
%         bw = imclearborder(bw, 8);
%         bw = bwmorph(bw, 'clean', Inf);
%         bw = bwmorph(bw, 'hbreak', Inf);
%         bw = bwmorph(bw, 'spur', Inf);
%         bw = imclose(bw, strel('disk', 4));
%         bw = imerode(bw, strel('disk', 2));
%         bw = bwmorph(bw, 'bridge', Inf);
%         bw = imdilate(bw, strel('disk', 3));
%         bw = bwpropfilt(bw, 'Eccentricity', [0 0.975]);
%         bw = imopen(bw, strel('disk', 20));
%         bw = imfill(bw, 'holes');
%     end
%     figure(); imshow(bw, [])

%     dosevalues = reshape((struct{i}.img_dose .* bw), 1, []);
%     dosevalues(dosevalues == 0) = [];
%     SEM = std(dosevalues(:))/sqrt(length(dosevalues)); % Standard Error
%     ts = tinv([0.025 0.975], length(dosevalues)-1); % T-Score
%     struct{i}.CI   = round(mean(dosevalues(:)) + ts*SEM, 1);
    struct{i}.mean = round(mean(struct{i}.ROI(:)), 1);
    struct{i}.std  = round(std(struct{i}.ROI(:)), 1);
%     figure(); imshow((struct{i}.img_dose .* bw), [])
    
    % Compute 95% CI for estimated dose
%     dosevalues = reshape(struct{i}.img_dose, 1, []);
%     dosevalues(dosevalues < 1.0) = []; 
%     dosevalues(dosevalues > 5.0) = [];
%     [min(dosevalues(:)) max(dosevalues(:))] 
%     
%     struct{i}.img_dose(struct{i}.img_dose > 1.0 & struct{i}.img_dose < 5.5) = 1;
% %     struct{i}.img_dose(struct{i}.img_dose <= 5.5) = 1;
%     struct{i}.img_dose(struct{i}.img_dose <= 1.0) = 0;
%     struct{i}.img_dose(struct{i}.img_dose >= 5.5) = 0;
%     
%     figure(); imshow(struct{i}.img_dose, [])
    
    % Plot    
    h = figure();
    hold on
    imshow(struct{i}.img_dose, [0 7], 'colormap', jet(4096))
    shading interp
    c = colorbar;
    c.Label.String = 'Dose (Gy)';
    title(strrep(struct{i}.filename, '_', '-'))
%     leg = legend(['D_{max}=' num2str()]);
%     title(leg, struct{i}.filename)
    set(gca, 'FontSize', 14)
%     saveas(h, sprintf('%s.png', fn))

    rectangle('Position', [ROI_bb{i} 100 100], 'EdgeColor','r', 'LineWidth', 3)
    
end

end