function [struct] = getEBT3struct(path, dpi, bit, ROI_size_mm, x_range, y_range)

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
channel = {'red', 'green', 'blue', 'gray'};

for i = 1:length(filelist)
    
    % Get directory, filename and extension of the image file
    [~, fn, e] = fileparts(filelist{i});
    
    struct{i}.filename = [fn e];
    
    % Read EBT3 film images
    [I_rgb, I_red, I_green, I_blue, I_gray] = readEBT3film(filelist{i});
    
    struct{i}.rgb.img   = I_rgb;
    struct{i}.red.img   = I_red;
    struct{i}.green.img = I_green;
    struct{i}.blue.img  = I_blue;
    struct{i}.gray.img  = I_gray;
    
    % Spatially register images to the first (control) image
    if i > 1
        tform = registerImage(struct{1}.gray.img, struct{i}.gray.img);
    end
    
    for j = 1:length(channel)
        
        if i > 1
            struct{i}.(channel{j}).img = imwarp( ...
                struct{i}.(channel{j}).img, tform, 'OutputView', ...
                imref2d(size(struct{1}.(channel{j}).img)));
        end
        
        % Crop and store the image segment which contains the selected ROI
        struct{i}.(channel{j}).ROI = imcrop(struct{i}.(channel{j}).img, ...
            ROI_bb);
        
        % Compute and store the mean and standard deviation of each ROI
        struct{i}.(channel{j}).PV    = mean2(struct{i}.(channel{j}).ROI);
        struct{i}.(channel{j}).sigma = std2(struct{i}.(channel{j}).ROI);
        
        % Compute optical density (OD)
        struct{i}.(channel{j}).OD = log10(2^bit/struct{i}.(channel{j}).PV);
        
    end
    
    % Plot superimposed ROI on EBT3 film
    figure(); imshow(struct{i}.gray.img, [])
    hold on
    rectangle('Position', ROI_bb, 'EdgeColor','r', 'LineWidth', 3)
     
    % Plot line profile (roughly) through the EBT3 film 
    x = [0 size(struct{i}.rgb.img,2)];
    y = [size(struct{i}.rgb.img,1)/2 size(struct{i}.rgb.img,1)/2];
    c = improfile(struct{i}.rgb.img, x, y);
    
    figure()
    subplot(2,1,1)
    imshow(struct{i}.rgb.img)
    hold on
    plot(x, y, 'r')
    subplot(2,1,2)
    plot(c(:,1,1), 'r')
    hold on
    plot(c(:,1,2), 'g')
    plot(c(:,1,3), 'b')
    hold off
    
end

end