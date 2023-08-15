function [struct] = getEBT3struct(path, dpi, bit, ROI_size_mm, x_range, y_range, flag)

%% Find and store all image files in directory
filelist = getAllFiles(path);

%% Compute ROI
px_size = 25.4/dpi;     % in mm/pixel
ROI_size_px = round(ROI_size_mm/px_size);

% Define center of ROI
x_center = (x_range(2) - x_range(1)) / 2;
y_center = ((y_range(2) - y_range(1)) / 2) + y_range(1);

% Coordinate of upper left corner of ROI
x_bb = round(x_center - ROI_size_px(1)/2);
y_bb = round(y_center - ROI_size_px(2)/2);

% Define ROI bounding box
ROI_bb = [x_bb y_bb ROI_size_px(1)-1 ROI_size_px(2)-1];

%% Loop through all files

% Initialization of image structures and vectors
struct = {};
channel = {'red', 'green', 'blue', 'gray', 'rgb'};

% [95 281] <-- ref!
% translation_vec = {[82 298]-[95 281], [95 281]-[93 259], ...
%     [95 281]-[100 262], [95 281]-[93 265], [95 281]-[82 264]};
translation_vec = {[-13 17], [-4 7], [-4 15], [22 2], [19 -5], [16 2], [17 13]};
counter = 0;

for i = 1:length(filelist)

    %     destpath = 'C:\Users\delmo\Desktop\Results\EBT3';
    if flag
        if ismember(i, [6 8 11:15])
            counter = counter + 1;
        end
    end

    % Get directory, filename and extension of the image file
    [p, fn, e] = fileparts(filelist{i});
    [~, p1, ~] = fileparts(p);

    struct{i}.filename = [fn e];

    % Read EBT3 film images
    [I_rgb, I_red, I_green, I_blue, I_gray] = readEBT3film(filelist{i});

    struct{i}.rgb.img   = I_rgb;

    if flag
        if ismember(i, 4:11)
            struct{i}.red.img   = flip(I_red,2);
            struct{i}.green.img = flip(I_green,2);
            struct{i}.blue.img  = flip(I_blue,2);
            struct{i}.gray.img  = flip(I_gray,2);
        else
            struct{i}.red.img   = I_red;
            struct{i}.green.img = I_green;
            struct{i}.blue.img  = I_blue;
            struct{i}.gray.img  = I_gray;
        end
    else
        struct{i}.red.img   = I_red;
        struct{i}.green.img = I_green;
        struct{i}.blue.img  = I_blue;
        struct{i}.gray.img  = I_gray;
    end

    % Spatially register images to the first (control) image
    if i > 1 % ~isequal(i,4) 

%         tform = imregcorr(struct{i}.gray.img, struct{1}.gray.img);

        [optimizer, metric] = imregconfig('multimodal');
        if flag
            optimizer.InitialRadius = 0.009;
            optimizer.Epsilon = 1.5e-4;
            optimizer.GrowthFactor = 1.01;
            optimizer.MaximumIterations = 500;
        end
        tform = imregtform(struct{i}.gray.img, struct{1}.gray.img, ...
            'rigid', optimizer, metric);

    end

    for j = 1:length(channel)

        if i > 1

%             struct{i}.(channel{j}).img = imwarp( ...
%                 struct{i}.(channel{j}).img, tform, 'OutputView', ...
%                 imref2d(size(struct{1}.(channel{j}).img)) );
%             [optimizer, metric] = imregconfig('multimodal');
%             optimizer.InitialRadius = 0.009;
%             optimizer.Epsilon = 1.5e-4;
%             optimizer.GrowthFactor = 1.01;
%             optimizer.MaximumIterations = 500;
%             struct{i}.(channel{j}).img = imregister( ...
%                 struct{i}.(channel{j}).img, struct{1}.(channel{j}).img, ...
%                 "affine", optimizer, metric, InitialTransformation = tform);

            struct{i}.(channel{j}).img = imwarp( ...
                struct{i}.(channel{j}).img, tform, 'OutputView', ...
                imref2d(size(struct{1}.(channel{j}).img)));

            if flag
                if ismember(i, [6 8 11:15])
                    struct{i}.(channel{j}).img = imtranslate( ...
                        struct{i}.(channel{j}).img, ...
                        translation_vec{counter});
                end
            end

        end

        if ~strcmp(channel{j}, 'rgb')

            % Median filtering
            %         struct{i}.(channel{j}).img = medfilt2( ...
            %             struct{i}.(channel{j}).img, [5 5]);

            % Crop and store the image segment which contains the selected
            % region of interest (ROI)
            struct{i}.(channel{j}).ROI = imcrop(struct{i}.(channel{j}).img, ...
                ROI_bb);

            % Compute and store the mean and standard deviation of ROI
            struct{i}.(channel{j}).PV    = mean2(struct{i}.(channel{j}).ROI);
            struct{i}.(channel{j}).sigma = std2(struct{i}.(channel{j}).ROI);

            % Compute optical density (OD) of ROI
            struct{i}.(channel{j}).OD = ...
                log10(2^bit./struct{i}.(channel{j}).PV);

            % Compute OD of the entire image
            struct{i}.(channel{j}).OD_img = ...
                log10(2^bit./struct{i}.(channel{j}).img);

        end

    end

    % Plot superimposed ROI on EBT3 film
%     figure(); imshow(struct{i}.gray.img, [])
%     hold on
%     rectangle('Position', ROI_bb, 'EdgeColor', 'r', 'LineWidth', 3)
%     title([strrep(struct{i}.filename, '_', '-') ...
%         ': grayscale, ' num2str(dpi) ' dpi'])
%     hold off

    %     destpath = fullfile(destpath, p1);
    %     mkdir(destpath);
    %     imwrite(uint16(struct{i}.rgb.img), fullfile(destpath, struct{i}.filename), ...
    %         'Compression', 'none', 'Resolution', 300)


    %     % Plot line profile (roughly) through the EBT3 film
    %     x = [0 size(struct{i}.rgb.img,2)];
    %     y = [size(struct{i}.rgb.img,1)/2 size(struct{i}.rgb.img,1)/2];
    %     struct{i}.x = [size(struct{i}.rgb.img,2)/2 size(struct{i}.rgb.img,2)/2];
    %     struct{i}.y = [0 size(struct{i}.rgb.img,1)];
    %     struct{i}.c = improfile(struct{i}.rgb.img, struct{i}.x, struct{i}.y);


    %     figure()
    %     subplot(2,1,1)
    %     imshow(struct{i}.rgb.img)
    %     hold on
    %     plot(x, y, 'r')
    %     subplot(2,1,2)
    %     plot(c(:,1,1), 'r')
    %     hold on
    %     plot(c(:,1,2), 'g')
    %     plot(c(:,1,3), 'b')
    %     hold off
    %     title(strrep(struct{i}.filename, '_', '-'))
    %     xlabel('y')
    %     ylabel('\it{PV}')
    %     set(gca, 'FontSize', 14)

end

end