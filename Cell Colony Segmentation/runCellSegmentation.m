clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%% Source and destination directory of the CFU images
% sourcepath_img_colonies = {'C:\Users\delmo\Desktop\OneDrive_1_8-30-2020'};
% sourcepath_img_colonies = { ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\T47D\X-ray\04122019', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\T47D\X-ray\18122019'};
% sourcepath_img_colonies = { ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\02122020', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\03012020', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\16122020', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\17122020', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\18112019', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\20112019', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\06022020', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\13022020', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\30012020'};
% sourcepath_img_colonies = { ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\Exp 1', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\Exp 7'};
% sourcepath_img_colonies = { ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\Proton\Tuesday2019', ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\Proton\Thursday2019'};
% sourcepath_img_colonies = {'C:\Users\delmo\Desktop\Raw Images'};

% ICPR2020
% sourcepath_img_colonies = {'C:\Users\delmo\Desktop\Benchmark\Images'}; 
% sourcepath_img_dish = {'C:\Users\delmo\Desktop\Benchmark\Masks'}; 

% Cell dish border directory
% % sourcepath_img_dish = { ...
% %     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB.tif', ...
% %     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyIH.tif', ...
% %     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyHS.tif', ...
% %     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyTS.tif', ...
% %     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB2.jpg', ...
% %     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB3.tif'};

% sourcepath_ctrl     = 'C:\Users\Delmon Arous\Desktop\Images\Dataset3\Control';
% sourcepath_treat    = 'C:\Users\Delmon Arous\Desktop\Images\Dataset3\Treatment';
% sourcepath_dish     = 'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\Binary Masks\BW.txt';
% destpath            = 'C:\Users\Delmon Arous\Desktop\Results\Dataset3';

% BiGART 2021
% sourcepath_treat    = {'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\17122020', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\18112019', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\A549\X-ray\20112019'};
sourcepath_treat = {'C:\Users\Delmon Arous\Desktop\Jacob\Test_v4\Test Images', ...
    'C:\Users\Delmon Arous\Desktop\Images\A549\X-ray\17122020', ...
    'C:\Users\Delmon Arous\Desktop\Images\A549\X-ray\18112019', ...
    'C:\Users\Delmon Arous\Desktop\Images\A549\X-ray\20112019', ...
    'C:\Users\Delmon Arous\Desktop\Images\03012020'};
destpath ='C:\Users\Delmon Arous\Desktop\Results\Test_Images_v6';

% sourcepath_dish     = 'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\Binary Masks\BW_v2.txt';
sourcepath_dish     = 'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB.tif';

%% Define parameters

% Pixel size in mm/pixel (inch in mm per dpi)
dpi = 1200; % 1200 for celler og 314 for bakterier
param.px_size = 25.4/dpi; % mm/px

% Magnus A549: [3 3] og 6
% Ingunn T47D: [3 3] og 6
% Hilde A549:  [5 5] og 30

% [80 8000]
% T47D: [4 4] og 6 (PCA), [5 5] og 10 watershed
%     param.gaussfilt_size_pca    = [2 2]; % gaussfilt_size_kmeans{i};
%     param.obrcbr_size_pca       = 25; % obrcbr_size_kmeans{i};
%     param.gaussfilt_size        = [3 3]; % gaussfilt_size_watershed{i};
%     param.obrcbr_size           = 6; % obrcbr_size_watershed{i};

% [800 35000]
% Ecoli control 2: [3 3] og 6, nyeste [2 2] og 5; extractBLOBs [3 3] og 90; PCA2
% Ecoli control 3: [5 5] og 8, nyeste [14 14] og 30; extractBLOBs [3 3] og 90; PCA1
% Ecoli control 9: [5 5] og 8, nyeste [10 10] og 30; extractBLOBs [3 3] og 90; PCA1

% [700 20000]
% Klebsiella pneumoniae control 1:  [5 5] og 10, nyeste [6 6] og 25; extractBLOBs [3 3] og 63; PCA1
% Klebsiella pneumoniae control 3:  [4 4] og 6, nyeste [2 2] og 20; extractBLOBs [3 3] og 30; PCA1
% Klebsiella pneumoniae light 2:    [5 5] og 20, nyeste [9 9] og 25; extractBLOBs [3 3] og 63; PCA1

% [2500 20000]
% Pseudomonas Aeruginosa control 1:  [4 4] og 6, nyeste [8 8] og 25; extractBLOBs [3 3] og 80; PCA1
% Pseudomonas Aeruginosa control 2:  [4 4] og 6, nyeste [8 8] og 25; extractBLOBs [3 3] og 80; PCA1
% Pseudomonas Aeruginosa light 1:    [4 4] og 6, nyeste [8 8] og 25; extractBLOBs [3 3] og 80; PCA1

% [500 5000]
% Steph Au control 2:   [3 3] og 6, nyeste [6 6] og 20; extractBLOBs [3 3] og 30; PCA2
% Steph Au light 1:     [3 3] og 6, nyeste [6 6] og 20; extractBLOBs [3 3] og 30; PCA2
% Steph Au light 3:     [3 3] og 6, nyeste [6 6] og 20; extractBLOBs [3 3] og 30; PCA2

% Magnus A549: area(1) = 140, area(2) = 8000
% Ingunn T47D: area(1) = 80,  area(2) = 8000
% Hilde A549:  area(1) = 500, area(2) = 8000

% Define a search interval for colony area
area = [70 5000]; % cancer cell mean diameter 20-30 um
param.area = [area(1)/2 area(1) area(2) area(2)*3];

% Gaussian filter size
param.gaussfilt_size_pca = [2 2]; % [2 2] % [4 4] egentlig 
param.gaussfilt_size = [5 5];
%gaussfilt_size = {[3 3], [3 3], [3 3]};

% Opening-Closing-by-Reconstruction structuring element size
param.obrcbr_size_pca = 40;  % 2 or 8

% Watershed segmentation vector
param.ws_vec = 0.15:0.01:0.45; % 0.25:0.01:0.37 hilde og bakterier!

% Magnus A549: 0.2 / 0.25
% Ingunn T47D: 0.2 / 0.15
% Hilde A549: 0.15 / 0.05 for < 10 Gy og 0.5 / 0.1 for 10 Gy

% Post-segmentation intensity threshold
param.intensity_thresh = 0.15;

% Cell culture container extractor
param.extractcontainer = 1;

% Intensive segmentation
param.intensiveseg = 1;

% Color image binary flag
param.rgb = 1;

% Control image binary flag
param.ctrl = exist('sourcepath_ctrl', 'var');

% Treatment image binary flag
% param.treat = 0; % for transfer learning
param.treat = exist('sourcepath_treat', 'var');

f = uifigure;

%% Loop through source directories

% if param.ctrl == 1
%     
%     eigvec_ctrl = initiateCellSegmentation(sourcepath_ctrl, ...
%         sourcepath_dish, destpath, param, f);
% 
% end

eigvec_ctrl = [0.6618 -0.5665 -0.4908; 0.5881 -0.0136 0.8087; 0.4648 ... 
    0.8238 -0.3241];

% param.treat = exist('sourcepath_treat', 'var'); % for transfer learning
% param.ctrl  = 0; 

for i = 4 % 2:4 % :length(sourcepath_treat)

    initiateCellSegmentation(sourcepath_treat{i}, sourcepath_dish, ...
        destpath, param, f);

%     eigvec_treat = initiateCellSegmentation(sourcepath_treat{i}, ...
%         sourcepath_dish, destpath, param, f, eigvec_ctrl);
    
%     if param.ctrl == 1 && param.treat == 1
%         
%         eigvec_treat = initiateCellSegmentation(sourcepath_treat{i}, ...
%             sourcepath_dish, destpath, param, f, eigvec_ctrl, sourcepath_ref);
%         
%     elseif param.ctrl == 0 && param.treat == 1
%         
%         eigvec_treat = initiateCellSegmentation(sourcepath_treat{i}, ...
%             sourcepath_dish, destpath, param, f, sourcepath_ref);
%         
%     end

end


%%

% if param.rgb == 1
%     
%     for i = 1:length(sourcepath_img_colonies)
%         
%         folderList = getAllFolders(sourcepath_img_colonies{i});
%         %folderList_dish = getAllFolders(sourcepath_img_dish{i});
%         
%         for j = 1:length(folderList)
%             [path1, foldername1, ~] = fileparts(folderList{j});
%             [path2, foldername2, ~] = fileparts(path1);
%             destpath = fullfile(destpath_img_colonies, foldername2);
%             if ~exist(destpath, 'dir')
%                 mkdir(destpath)
%             end
%             destpath = fullfile(destpath, foldername1);
%             if ~exist(destpath, 'dir')
%                 mkdir(destpath)
%             end
%             
%             fileList = getAllFiles(folderList{j});
%             %fileList_dish = getAllFiles(folderList_dish{j});
%             
%             filename_vec = {};
%             confluency_vec = [];
%             n_colony = [];
%             mean_area   = []; SD_area   = [];
%             mean_red    = []; SD_red    = [];
%             mean_green  = []; SD_green  = [];
%             mean_blue   = []; SD_blue   = [];
%             mean_gray   = []; SD_gray   = [];
%             mean_pca1   = []; SD_pca1   = [];
%             mean_pca2   = []; SD_pca2   = [];
%             
%             for k = 1:length(fileList) % 2!!!!
%                 
%                 [~, filename, ~] = fileparts(fileList{k});
%                 filename
%                 
% %                 param.gaussfilt_size = gaussfilt_size{k};
% %                 param.obrcbr_size = obrcbr_size(k);
%                 
%                 t = tic;
%                 [confluency, areastats, redstats, greenstats, bluestats, ...
%                     graystats, pca1stats, pca2stats, ~, eigvec] = ...
%                     CellSegmentation(fileList{k}, sourcepath_img_dish{1}, ...
%                     destpath, param, f); % fileList_dish{k}
%                 runtime = toc(t);
%                 
%                 if ~isempty(areastats)
%                     
%                     colony_id           = [1:numel(areastats)]'; 
%                     colony_xyCentroids  = vertcat(areastats.Centroid);
%                     colony_xCentroids   = colony_xyCentroids(:,1);
%                     colony_yCentroids   = colony_xyCentroids(:,2);
%                     colony_areas        = [areastats.Area]' .* ...
%                         param.px_size^2;
%                     colony_circularity  = [areastats.Circularity]';
%                     colony_eccentricity = [areastats.Eccentricity]';
%                     confluency_vec(k)   = confluency * 100;
%                     
%                 else
%                     
%                     colony_id           = 0;
%                     colony_xCentroids   = 0;
%                     colony_yCentroids   = 0;
%                     colony_areas        = 0;
%                     colony_circularity  = 0;
%                     colony_eccentricity = 0;
%                     confluency_vec(k)   = 0;
%                     
%                 end
%                 
%                 colony_MeanRed      = [redstats.MeanIntensity]';
%                 colony_MeanGreen    = [greenstats.MeanIntensity]';
%                 colony_MeanBlue     = [bluestats.MeanIntensity]';
%                 colony_MeanGray     = [graystats.MeanIntensity]';
%                 colony_MeanPCA1     = [pca1stats.MeanIntensity]';
%                 colony_MeanPCA2     = [pca2stats.MeanIntensity]';
%                 colony_SDRed        = [redstats.SD]';
%                 colony_SDGreen      = [greenstats.SD]';
%                 colony_SDBlue       = [bluestats.SD]';
%                 colony_SDGray       = [graystats.SD]';
%                 colony_SDPCA1       = [pca1stats.SD]';
%                 colony_SDPCA2       = [pca2stats.SD]';
%                 
%                 filename_vec{k}         = filename;
%                 n_colony(k)             = numel(areastats);
%                 mean_area(k)            = mean(colony_areas);
%                 mean_circularity(k)     = mean(colony_circularity);
%                 mean_eccentricity(k)    = mean(colony_eccentricity);
%                 mean_red(k)             = mean(colony_MeanRed);
%                 mean_green(k)           = mean(colony_MeanGreen);
%                 mean_blue(k)            = mean(colony_MeanBlue);
%                 mean_gray(k)            = mean(colony_MeanGray);
%                 mean_pca1(k)            = mean(colony_MeanPCA1);
%                 mean_pca2(k)            = mean(colony_MeanPCA2);
%                 SD_area(k)              = std(colony_areas);
%                 SD_circularity(k)       = std(colony_circularity);
%                 SD_eccentricity(k)      = std(colony_eccentricity);
%                 SD_red(k)               = mean(colony_SDRed);
%                 SD_green(k)             = mean(colony_SDGreen);
%                 SD_blue(k)              = mean(colony_SDBlue);
%                 SD_gray(k)              = mean(colony_SDGray);
%                 SD_pca1(k)              = mean(colony_SDPCA1);
%                 SD_pca2(k)              = mean(colony_SDPCA2);
%                 
%                 % pixelindexlist{k} = cat(1,areastats_colony.PixelIdxList);
%                 
%                 % Write results to .csv file
%                 varNames = {'Colony ID',            ...
%                     'Colony Area (mm2)',            ...
%                     'Centroid x-Coordinate (px)',   ...
%                     'Centroid y-Coordinate (px)',   ...
%                     'Circularity (a.u)',            ...
%                     'Eccentricity (a.u)',           ...
%                     'Mean Colony Red (a.u)',        ...
%                     'Mean Colony Green (a.u)',      ...
%                     'Mean Colony Blue (a.u)',       ...
%                     'Mean Colony Gray (a.u)',       ...
%                     'Mean Colony PCA1 (a.u)',       ...
%                     'Mean Colony PCA2 (a.u)',       ...
%                     'SD Colony Red (a.u)',          ...
%                     'SD Colony Green (a.u)',        ...
%                     'SD Colony Blue (a.u)',         ...
%                     'SD Colony Gray (a.u)',         ...
%                     'SD Colony PCA1 (a.u)',         ...
%                     'SD Colony PCA2 (a.u)'};
%                 T = table(colony_id,        ...
%                     colony_areas,           ...
%                     colony_xCentroids,      ...
%                     colony_yCentroids,      ...
%                     colony_circularity,     ...
%                     colony_eccentricity,    ...
%                     colony_MeanRed,         ...
%                     colony_MeanGreen,       ...
%                     colony_MeanBlue,        ...
%                     colony_MeanGray,        ...
%                     colony_MeanPCA1,        ...
%                     colony_MeanPCA2,        ...
%                     colony_SDRed,           ...
%                     colony_SDGreen,         ...
%                     colony_SDBlue,          ...
%                     colony_SDGray,          ...
%                     colony_SDPCA1,          ...
%                     colony_SDPCA2,          ...
%                     'VariableNames', varNames);
%                 writetable(T, fullfile(destpath, ...
%                     sprintf('%s-ColonyData.csv', filename)),  ...
%                     'Delimiter', ',', 'QuoteStrings', true)
%                 
%             end
%             
%             varNames  = {'Filename',        ...
%                 'Colony count',             ...
%                 'Mean area (mm2)',      	...
%                 'SD area (mm2)',            ...
%                 'Confluency (%)',           ...
%                 'Mean circularity (a.u)',   ...
%                 'SD circularity (a.u)',     ...
%                 'Mean eccentricity (a.u)',  ...
%                 'SD eccentricity (a.u)',    ...
%                 'Mean Red (a.u)',           ...
%                 'SD Red (a.u)',             ...
%                 'Mean Green (a.u)',         ...
%                 'SD Green (a.u)',           ...
%                 'Mean Blue (a.u)',          ...
%                 'SD Blue (a.u)',            ...
%                 'Mean Gray (a.u)',          ...
%                 'SD Gray (a.u)',            ...
%                 'Mean PCA1 (a.u)',          ...
%                 'SD PCA1 (a.u)',            ...
%                 'Mean PCA2 (a.u)',          ...
%                 'SD PCA2 (a.u)'};
%             T = table(filename_vec',    ...
%                 n_colony',              ...
%                 mean_area',             ...
%                 SD_area',               ...
%                 confluency_vec',        ...
%                 mean_circularity',      ...
%                 SD_circularity',        ...
%                 mean_eccentricity',     ...
%                 SD_eccentricity',       ...
%                 mean_red',              ...
%                 SD_red',                ...
%                 mean_green',            ...
%                 SD_green',              ...
%                 mean_blue',             ...
%                 SD_blue',               ...
%                 mean_gray',             ...
%                 SD_gray',               ...
%                 mean_pca1',             ...
%                 SD_pca1',               ...
%                 mean_pca2',             ...
%                 SD_pca2',               ...
%                 'VariableNames', varNames);
%             writetable(T, fullfile(destpath, [foldername2 '.csv']),  ...
%                 'Delimiter', ',', 'QuoteStrings', true)
%             
%             close all
%             
%         end
%         
%     end
%     
% else
%     
%     for i = 1:length(sourcepath_img_colonies)
%         
%         folderList = getAllFolders(sourcepath_img_colonies{i});
% %         folderList_dish = getAllFolders(sourcepath_img_dish{i});
%         
%         for j = 1:length(folderList)
%             [path1, foldername1, ~] = fileparts(folderList{j});
%             [path2, foldername2, ~] = fileparts(path1);
%             destpath = fullfile(destpath_img_colonies, foldername2);
%             if ~exist(destpath, 'dir')
%                 mkdir(destpath)
%             end
%             destpath = fullfile(destpath, foldername1);
%             if ~exist(destpath, 'dir')
%                 mkdir(destpath)
%             end
%             
%             fileList = getAllFiles(folderList{j});
% %             fileList_dish = getAllFiles(folderList_dish{j});
%             
%             filename_vec = {};
%             confluency_vec = [];
%             n_colony = [];
%             mean_area       = []; SD_area       = [];
%             mean_intensity  = []; SD_intensity  = [];
%             mean_pca1       = []; SD_pca1       = [];
%             
%             for k = 1:length(fileList)
%                 
%                 [~, filename, ~] = fileparts(fileList{k});
%                 filename
%                 
%                 t = tic;
%                 [confluency, areastats, imgstats, pca1stats] = ...
%                     CellSegmentation_test(fileList{k}, ...
%                     sourcepath_img_dish{6}, destpath, param, f); % fileList_dish{k} 
%                 runtime = toc(t);
%                 
%                 if ~isempty(areastats)
%                     
%                     colony_xyCentroids = vertcat(areastats.Centroid);
%                     colony_xCentroids = colony_xyCentroids(:,1);
%                     colony_yCentroids = colony_xyCentroids(:,2);
%                     colony_areas = [areastats.Area]' .* param.px_size^2;
%                     confluency_vec(k) = confluency * 100;
%                     
%                 else
%                     
%                     colony_xCentroids = 0;
%                     colony_yCentroids = 0;
%                     colony_areas = 0;
%                     confluency_vec(k) = 0;
%                     
%                 end
%                 
%                 colony_MedianIntensity  = [imgstats.MedianIntensity]';
%                 colony_MedianPCA1       = [pca1stats.MedianIntensity]';
%                 colony_SDIntensity      = [imgstats.SD]';
%                 colony_SDPCA1           = [pca1stats.SD]';
%                 
%                 filename_vec{k}     = filename;
%                 n_colony(k)         = numel(areastats);
%                 mean_area(k)        = mean(colony_areas);
%                 mean_intensity(k)   = mean(colony_MedianIntensity);
%                 mean_pca1(k)        = mean(colony_MedianPCA1);
%                 SD_area(k)          = mean(colony_areas);
%                 SD_intensity(k)     = mean(colony_SDIntensity);
%                 SD_pca1(k)          = mean(colony_SDPCA1);
%                 
%                 % pixelindexlist{k} = cat(1,areastats_colony.PixelIdxList);
%                 
%                 % Write results to .csv file
%                 varNames = {'Colony Area (mm2)',        ...
%                     'Centroid x-Coordinate (px)',       ...
%                     'Centroid y-Coordinate (px)',       ...
%                     'Median Colony Intensity (a.u)',    ...
%                     'Median Colony PCA1 (a.u)',         ...
%                     'SD Colony Intensity (a.u)',        ...
%                     'SD Colony PCA1 (a.u)'};
%                 T = table(colony_areas, ...
%                     colony_xCentroids,  ...
%                     colony_yCentroids,  ...
%                     colony_MedianIntensity, ...
%                     colony_MedianPCA1,  ...
%                     colony_SDIntensity,     ...
%                     colony_SDPCA1,      ...
%                     'VariableNames', varNames);
%                 writetable(T, fullfile(destpath, ...
%                     sprintf('%s-ColonyData.csv', filename)),  ...
%                     'Delimiter', ',', 'QuoteStrings', true)
%                 
%             end
%             
%             varNames  = {'Filename',    ...
%                 'Colony count',         ...
%                 'Mean size (mm2)',      ...
%                 'SD size (mm2)',        ...
%                 'Confluency (%)',       ...
%                 'Mean Intensity (a.u)', ...
%                 'SD Intensity (a.u)',   ...
%                 'Mean PCA1 (a.u)',      ...
%                 'SD PCA1 (a.u)'};
%             T = table(filename_vec',    ...
%                 n_colony',              ...
%                 mean_area',             ...
%                 SD_area',               ...
%                 confluency_vec',        ...
%                 mean_intensity',        ...
%                 SD_image',              ...
%                 mean_pca1',             ...
%                 SD_pca1',               ...
%                 'VariableNames', varNames);
%             writetable(T, fullfile(destpath, [foldername2 '.csv']),  ...
%                 'Delimiter', ',', 'QuoteStrings', true)
%             
%             close all
%             
%         end
%         
%     end
%     
% end

% delete(gcp)
