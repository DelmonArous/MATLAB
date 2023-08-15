clear all
close all
clc

%%
pathCT  = 'C:\Users\delmo\Desktop\DICOMreslice_v8_test - Copy';
path_bw = 'C:\Users\delmo\Desktop\Mask.txt';

%% Loop through all DICOM CT images
filelist    = getAllFiles(pathCT);
filenames   = {};

for i = 1:length(filelist)
    
    [path, name, ext] = fileparts(filelist{i});
    filenames{i} = [name ext];
    
end

filenames = filenames.';

%% Extract DICOM CT images
imagesCT = readDICOMImages(pathCT, filenames);

size(imagesCT.data)

% Re-order the slices
% ENTEN FLIP(X,1) eller FLIP(X,3)
imagesCT.data = flip(imagesCT.data, 3);

size(imagesCT.data)

% Set image widths and start values
width = [imagesCT.width(1) imagesCT.width(2)];
start = [imagesCT.start(1) imagesCT.start(2)];

%% Read binary mask
img_bw  = load(path_bw);
bw      = double(~img_bw);
bw      = (bw - min(min(bw))) ./ abs(max(max(bw)) - min(min(bw)));

figure(); imshow(bw, [])

%% Loop over CT slices

for slice = 36:55
    
    imgCT_temp  = imagesCT.data(:,:,slice);
    
    % Remove the extra dimension
    imgCT_temp = squeeze(imgCT_temp)';
    
    % Pad the image and calculate the start offsets
    s = max(size(imgCT_temp') .* width);
    offset = ((size(imgCT_temp') .* width) - s) / 2;
    
    % Create spatial reference object based on the start and width inputs
    RI = imref2d(size(imgCT_temp), [start(1) + offset(1) start(1) + ...
        offset(1) + size(imgCT_temp,2) * width(1)], [-(start(2) - offset(2) + ...
        size(imgCT_temp,1) * width(2)) -start(2) + offset(2)]);
    
    % Plot
    h = figure();
%     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(ind2rgb(gray2ind((imgCT_temp) / 2048, 256), colormap('gray')), RI)
    xlabel('cm')
    ylabel('cm')
    set(gca, 'FontSize', 14)
    
    
end


