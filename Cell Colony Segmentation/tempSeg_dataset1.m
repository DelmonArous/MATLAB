clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%%
sourcepath  = {'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Manuscripts\Colony Counter paper\Images\Dataset1'};
% sourcepath = {'C:\Users\Delmon Arous\Desktop\Images\Dataset2\E coli', ...
%     'C:\Users\Delmon Arous\Desktop\Images\Dataset2\Klebsiella pneumoniae', ...
%     'C:\Users\Delmon Arous\Desktop\Images\Dataset2\Pseudomonas aeruginosa', ...
%     'C:\Users\Delmon Arous\Desktop\Images\Dataset2\Staph Au'};
destpath    = 'C:\Users\delmo\Desktop\Results\Dataset1';

% Cell dish border directory
% sourcepath_img_dish = { ...
%     'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB.tif', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyIH.tif', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyHS.tif', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyTS.tif', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB2.jpg', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\EmptyMB3.tif'};

% sourcepath_img_dish = { ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Colony Counter paper\Results ICPR2020\Dish Template Masks\E coli masks', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Colony Counter paper\Results ICPR2020\Dish Template Masks\Klebsiella pneumoniae masks', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Colony Counter paper\Results ICPR2020\Dish Template Masks\Pseudomonas aeruginosa masks', ...
%     'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Colony Counter paper\Results ICPR2020\Dish Template Masks\Staph Au masks'};
sourcepath_img_dish = {'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\Binary Masks\BW.txt'};

%% Define parameters

% Pixel size in mm/pixel (inch in mm per dpi)
dpi = 1200; % 1200 for celler og 314 for bakterier
param.px_size = 25.4/dpi;

% Define a search interval for colony area
area = [60 7000]; % 76 for T47D cancer cell mean diameter 20-30 um
param.area = [area(1)/2 area(1) area(2) area(2)*2];

% Gaussian filtering and Opening-Closing-by-Reconstruction prior to k-Means
% param.gaussfilt_size_pca = [3 3];
% param.obrcbr_size_pca = 30; 
gaussfilt_size_kmeans = {[4 4], [3 3], [4 4]};
obrcbr_size_kmeans = {55, 20, 55};

% Gaussian filtering and Opening-Closing-by-Reconstruction prior to
% watershed
% param.gaussfilt_size = [2 2];
% param.obrcbr_size = 4;
gaussfilt_size_watershed = {[5 5], [2 2], [5 5]};
obrcbr_size_watershed = {25, 20, 25};

% Top-hat filtering size
param.tophat_size = 90;

% Watershed segmentation vector
param.ws_vec = 0.15:0.01:0.37;

% Post-segmentation intensity threshold
param.intensity_thresh = 0.1; % 0.1 for T47D

% Cell culture container extractor
param.extractcontainer = 1;

% Intensive segmentation
param.intensiveseg = 0;

% Color image flag
param.rgb = 1;

%%

colors          = ['r' 'g' 'b']; %  , 'c', 'm', 'y'};
filename_vec    = {};

filelist        = getAllFiles(sourcepath{1});
% filelist_dish   = getAllFiles(sourcepath_img_dish{2});
f = uifigure;

for i = 1:length(filelist)
    
    [path, filename, ext] = fileparts(filelist{i});
    param.gaussfilt_size_pca    = [2 2]; % gaussfilt_size_kmeans{i};     %  
    param.obrcbr_size_pca       = 40;    % obrcbr_size_kmeans{i};        %       
    param.gaussfilt_size        = [3 3]; % gaussfilt_size_watershed{i};  %   
    param.obrcbr_size           = 6;     % obrcbr_size_watershed{i};     %    
    
    if usejava('jvm') && feature('ShowFigureWindows')
        
        % Start progress bar
        progress = uiprogressdlg(f, 'Title', ...
            'Segmenting colonies...', 'Message', '1', 'Cancelable', 'on');
        progress.Value = 0.0;
        progress.Message = sprintf( ...
            '%s%s:\nReading colony image', ...
            filename, ext);
        
    end
    
    % Attempt to load image file using imread
    try
        
        % If imread is successful, store the image information
        [img_in, ~] = imread(filelist{i});
        img_original = img_in(:,:,1:3);
        
    catch
        
        % Otherwise, the file is either corrupt or not a real image.
        % Throw an error and stop execution
        disp(['File ', filename, ' is not a valid cell colony image!']);
        
    end
    
    % Dimensions of input image data array
    [rows, cols, dims] = size(img_in);
    
    % Check if the image array contains RGB channels. If so, extract RGB color
    % channels and normalize them between [0,1]
    if dims == 3 || dims == 4
        
        % Store red, green and blue channel, respectively
        redchannel = double(img_in(:,:,1));
        greenchannel = double(img_in(:,:,2));
        bluechannel = double(img_in(:,:,3));
        img_in = img_in(:,:,1:3);
        
        % Convert the RGB (truecolor) image to grayscale
        img_in = double(rgb2gray(img_in));
        
        % Min-max normalization of each color channel image
        img_red = (redchannel - min(min(redchannel))) ./ ...
            abs(max(max(redchannel)) - min(min(redchannel)));
        img_green = (greenchannel - min(min(greenchannel))) ./ ...
            abs(max(max(greenchannel)) - min(min(greenchannel)));
        img_blue = (bluechannel - min(min(bluechannel))) ./ ...
            abs(max(max(bluechannel)) - min(min(bluechannel)));
        
    elseif dims < 3
        
        % Otherwise, use the first (red) channel as the reference image
        img_in = double(img_in(:,:,1));
        
    else
        
        % Throw an error and stop execution
        disp('The image data array contains too many dimensions.');
        
    end
    
    % Normalization
    [img_norm, img_norm2] = normalizationMinMax(img_in);
    
    % Store normalized image (either grayscale or red channel)
    img_plane = img_norm;
    
    % Mask out the inner space within the dish/flask
    if param.extractcontainer == 1
        bw_border = getCellContainerBorder(sourcepath_img_dish{1}, ... % filelist_dish{i}
            img_norm, f);
    else
        bw_border = ones(rows, cols);
    end
    
    %     figure(); imshow(img_norm, [])
%         figure(); imshow(bw_border, [])
    
    % Multiply the dish/flask mask with the colony image to extract out
    % solely free viable area in the container
    img_norm = bw_border .* img_norm;
    
%     figure(); imshow(img_norm, [])
    
    %% Independent component analysis (ICA)
    
    % Perform Fast ICA
%     [Z_fastica, W, T, mu] = fastICA( ...
%         double(reshape(img_original, rows * cols, dims)).', dims, ...
%         'kurtosis', 1);
%     Z_fastica = Z_fastica.';
%     fastICAstruct = ChannelSelection(Z_fastica, rows, cols, dims, ...
%         bw_border, 'fastICA', filename, destpath);
    
%     h_fastica1 = figure(); imshow(fastICAstruct.channel1.img, []); title('ICA1 fast-kurtosis')
%     h_fastica2 = figure(); imshow(fastICAstruct.channel2.img, []); title('ICA2 fast-kurtosis')
%     h_fastica3 = figure(); imshow(fastICAstruct.channel3.img, []); title('ICA3 fast-kurtosis')
%     saveas(h_fastica1, fullfile(destpath, ...
%         sprintf('%s_fastkurt-ICA1.fig', filename)))
%     saveas(h_fastica2, fullfile(destpath, ...
%         sprintf('%s_fastkurt-ICA2.fig', filename)))
%     saveas(h_fastica3, fullfile(destpath, ...
%         sprintf('%s_fastkurt-ICA3.fig', filename)))
    
    % Perform max-kurtosis ICA
%     [Z_maxica, W, T, mu] = kICA( ...
%         double(reshape(img_original, rows * cols, dims)).', dims);
%     Z_maxica = Z_maxica.';
%     maxICAstruct = ChannelSelection(Z_maxica, rows, cols, dims, ...
%         bw_border, 'maxICA', filename, destpath);
%     
%     h_maxica1 = figure(); imshow(maxICAstruct.channel1.img, []); title('ICA1 max-kurtosis')
%     h_maxica2 = figure(); imshow(maxICAstruct.channel2.img, []); title('ICA2 max-kurtosis')
%     h_maxica3 = figure(); imshow(maxICAstruct.channel3.img, []); title('ICA3 max-kurtosis')
%     saveas(h_maxica1, fullfile(destpath, ...
%         sprintf('%s_maxkurt-ICA1.fig', filename)))
%     saveas(h_maxica2, fullfile(destpath, ...
%         sprintf('%s_maxkurt-ICA2.fig', filename)))
%     saveas(h_maxica3, fullfile(destpath, ...
%         sprintf('%s_maxkurt-ICA3.fig', filename)))
    
    %% Histogram plot
    %     h_fasticahist = figure();
    %     hold on
    %     histogram(img_fastica1(:), 2*256);
    %     histogram(img_fastica2(:), 2*256);
    %     histogram(img_fastica3(:), 2*256);
    %     hold off
    %     legend('ICA1', 'ICA2', 'ICA3')
    %     grid on;
    %     xlim([-5 5])
    %     saveas(h_fasticahist, fullfile(destpath, ...
    %         sprintf('%s_ICAfastkurt-hist.fig', filename)))
    
    %     h_maxicahist = figure();
    %     hold on
    %     histogram(img_maxica1(:), 2*256);
    %     histogram(img_maxica2(:), 2*256);
    %     histogram(img_maxica3(:), 2*256);
    %     hold off
    %     legend('ICA1', 'ICA2', 'ICA3')
    %     grid on;
    %     xlim([-5 5])
    %     saveas(h_maxicahist, fullfile(destpath, ...
    %         sprintf('%s_ICAmaxkurt-hist.fig', filename)))
    
    %% Principal Component Analysis (PCA)
    
    if dims == 3 || dims == 4
        
        % Perform PCA on the three bundled color channels
        [~, Z_pca, ~, explained, img_pca1, img_pca2, img_pca3] = ...
            runPCA(double(img_original));
        
        pca_explained_thresh = 0.000000001;
        
        % The clonogenic information is usually projected onto the
        % 2nd principal component plane
        if ~isempty(img_pca2) && explained(2) > pca_explained_thresh
            
            % Therefore, store the PCA2 image as the primary PCA image
            img_pca = img_pca2;
            param.pca = 'PCA2';
            
        elseif (isempty(img_pca2) || explained(2) <= pca_info_thresh) ...
                && ~isempty(img_pca1)
            
            % Otherwise, not enough color variation is provided in the
            % input image such that no information is projected onto the
            % 2nd principal component plane, but rather onto the 1st
            % prinicpal component plane. Therefore, store the PCA1 image
            % as the PCA image
            img_pca = img_pca1;
            param.pca = 'PCA1';
            
        else
            
            % Throw an error and stop execution
            disp('Principal component analysis (PCA) failed!')
            
        end
        
        % Otherwise, the input image consist of one
    else
        
        % Perform PCA on the single image channel
        [~, Z_pca, ~, ~, img_pca1] = runPCA(double(img_original));
        
        % The clonogenic information is now projected onto the
        % 1st principal component plane
        if ~isempty(img_pca1)
            
            % Therefore, store the PCA1 image as the primary PCA image
            img_pca = img_pca1;
            param.pca = 'PCA1';
            
        else
            
            % Throw an error and stop execution
            disp('Principal component analysis (PCA) failed!')
            
        end
        
    end
    
    img_pca = img_pca1;
    PCAstruct = ChannelSelection(Z_pca, rows, cols, dims, ...
        bw_border, 'PCA', filename, destpath);
%     PCAstruct.channelopt.img = img_pca1;
    
%     h_PCA1 = figure(); imshow(PCAstruct.channel1.img, []); title('PCA1')
%     h_PCA2 = figure(); imshow(PCAstruct.channel2.img, []); title('PCA2')
%     h_PCA3 = figure(); imshow(PCAstruct.channel3.img, []); title('PCA3')
%     saveas(h_PCA1, fullfile(destpath, sprintf('%s_PCA1.fig', filename)))
%     saveas(h_PCA2, fullfile(destpath, sprintf('%s_PCA2.fig', filename)))
%     saveas(h_PCA3, fullfile(destpath, sprintf('%s_PCA3.fig', filename)))
    
    %%
    filename_vec{i}     = filename;
    
    % Store Fast Kurtosis ICA results
%     contrast_fastica1(i)     = fastICAstruct.channel1.Contrast;
%     homogeneity_fastica1(i)  = fastICAstruct.channel1.Homogeneity;
%     correlation_fastica1(i)  = fastICAstruct.channel1.Correlation;
%     energy_fastica1(i)       = fastICAstruct.channel1.Energy;
%     entropy_fastica1(i)      = fastICAstruct.channel1.Entropy;
%     SNR_fastica1(i)          = fastICAstruct.channel1.SNR;
%     Otsu_fastica1(i)         = fastICAstruct.channel1.OtsuSigma;
%     Log_fastica1(i)          = fastICAstruct.channel1.threshOutLog;
%     
%     contrast_fastica2(i)     = fastICAstruct.channel2.Contrast;
%     homogeneity_fastica2(i)  = fastICAstruct.channel2.Homogeneity;
%     correlation_fastica2(i)  = fastICAstruct.channel2.Correlation;
%     energy_fastica2(i)       = fastICAstruct.channel2.Energy;
%     entropy_fastica2(i)      = fastICAstruct.channel2.Entropy;
%     SNR_fastica2(i)          = fastICAstruct.channel2.SNR;
%     Otsu_fastica2(i)         = fastICAstruct.channel2.OtsuSigma;
%     Log_fastica2(i)          = fastICAstruct.channel2.threshOutLog;
%     
%     contrast_fastica3(i)     = fastICAstruct.channel3.Contrast;
%     homogeneity_fastica3(i)  = fastICAstruct.channel3.Homogeneity;
%     correlation_fastica3(i)  = fastICAstruct.channel3.Correlation;
%     energy_fastica3(i)       = fastICAstruct.channel3.Energy;
%     entropy_fastica3(i)      = fastICAstruct.channel3.Entropy;
%     SNR_fastica3(i)          = fastICAstruct.channel3.SNR;
%     Otsu_fastica3(i)         = fastICAstruct.channel3.OtsuSigma;
%     Log_fastica3(i)          = fastICAstruct.channel3.threshOutLog;
    
%     % Store Max Kurtosis ICA results
%     contrast_maxica1(i)     = maxICAstruct.channel1.GLCM.Contrast;
%     homogeneity_maxica1(i)  = maxICAstruct.channel1.GLCM.Homogeneity;
%     correlation_maxica1(i)  = maxICAstruct.channel1.GLCM.Correlation;
%     energy_maxica1(i)       = maxICAstruct.channel1.GLCM.Energy;
%     kurtosis_maxica1(i)     = maxICAstruct.channel1.kurtosis;
%     SNR_maxica1(i)          = maxICAstruct.channel1.SNR;
%     Otsu_maxica1(i)         = maxICAstruct.channel1.OtsuSigma;
%     Log_maxica1(i)          = maxICAstruct.channel1.threshOutLog;
%     
%     contrast_maxica2(i)     = maxICAstruct.channel2.GLCM.Contrast;
%     homogeneity_maxica2(i)  = maxICAstruct.channel2.GLCM.Homogeneity;
%     correlation_maxica2(i)  = maxICAstruct.channel2.GLCM.Correlation;
%     energy_maxica2(i)       = maxICAstruct.channel2.GLCM.Energy;
%     kurtosis_maxica2(i)     = maxICAstruct.channel2.kurtosis;
%     SNR_maxica2(i)          = maxICAstruct.channel2.SNR;
%     Otsu_maxica2(i)         = maxICAstruct.channel2.OtsuSigma;
%     Log_maxica2(i)          = maxICAstruct.channel2.threshOutLog;
%     
%     contrast_maxica3(i)     = maxICAstruct.channel3.GLCM.Contrast;
%     homogeneity_maxica3(i)  = maxICAstruct.channel3.GLCM.Homogeneity;
%     correlation_maxica3(i)  = maxICAstruct.channel3.GLCM.Correlation;
%     energy_maxica3(i)       = maxICAstruct.channel3.GLCM.Energy;
%     kurtosis_maxica3(i)     = maxICAstruct.channel3.kurtosis;
%     SNR_maxica3(i)          = maxICAstruct.channel3.SNR;
%     Otsu_maxica3(i)         = maxICAstruct.channel3.OtsuSigma;
%     Log_maxica3(i)          = maxICAstruct.channel3.threshOutLog;
    
    % Store PCA results
    contrast_pca1(i)        = PCAstruct.channel1.Contrast;
    homogeneity_pca1(i)     = PCAstruct.channel1.Homogeneity;
    correlation_pca1(i)     = PCAstruct.channel1.Correlation;
    energy_pca1(i)          = PCAstruct.channel1.Energy;
    rawentropy_pca1(i)      = PCAstruct.channel1.rawEntropy;
    meanentropy_pca1(i)     = PCAstruct.channel1.meanEntropy;
    medianentropy_pca1(i)   = PCAstruct.channel1.medianEntropy;
    explained_pca1(i)       = explained(1);
    
    contrast_pca2(i)        = PCAstruct.channel2.Contrast;
    homogeneity_pca2(i)     = PCAstruct.channel2.Homogeneity;
    correlation_pca2(i)     = PCAstruct.channel2.Correlation;
    energy_pca2(i)          = PCAstruct.channel2.Energy;
    rawentropy_pca2(i)      = PCAstruct.channel2.rawEntropy;
    meanentropy_pca2(i)     = PCAstruct.channel2.meanEntropy;
    medianentropy_pca2(i)   = PCAstruct.channel2.medianEntropy;
    explained_pca2(i)       = explained(2);
    
    contrast_pca3(i)        = PCAstruct.channel3.Contrast;
    homogeneity_pca3(i)     = PCAstruct.channel3.Homogeneity;
    correlation_pca3(i)     = PCAstruct.channel3.Correlation;
    energy_pca3(i)          = PCAstruct.channel3.Energy;
    rawentropy_pca3(i)      = PCAstruct.channel3.rawEntropy;
    meanentropy_pca3(i)     = PCAstruct.channel3.meanEntropy;
    medianentropy_pca3(i)   = PCAstruct.channel3.medianEntropy;
    explained_pca3(i)       = explained(3);
    
    %% Nå sett kriterium basert på metrics og velg én kanal fra hver
    
    PCAstruct.channelopt.img = img_norm;
   
%     bw_BLOBs_fastICA = extractBLOBs(fastICAstruct.channelopt.img, bw_border, param);
%     bw_BLOBs_maxICA = extractBLOBs(maxICAstruct.channelopt.img, bw_border, param);    
    bw_BLOBs_PCA = extractBLOBs(PCAstruct.channelopt.img, bw_border, param);
    
%     h_bw_BLOBs_fastICA = figure(); imshow(bw_BLOBs_fastICA, []); title('ICA fast-kurtosis')
%     h_bw_BLOBs_maxICA = figure(); imshow(bw_BLOBs_maxICA, []); title('ICA max-kurtosis')
    h_bw_BLOBs_PCA = figure(); imshow(bw_BLOBs_PCA, []); title('PCA')
%     saveas(h_bw_BLOBs_fastICA, fullfile(destpath, sprintf('%s_bwBLOBsFastICA.fig', filename)))
%     saveas(h_bw_BLOBs_maxICA, fullfile(destpath, sprintf('%s_bwBLOBsMaxICA.fig', filename)))
    saveas(h_bw_BLOBs_PCA, fullfile(destpath, sprintf('%s_bwBLOBsPCA.fig', filename)))
    
    %% Guassian smoothing and removal of outliers
    
    % Gaussian filtering
    img_norm1 = imgaussfilt(img_norm, param.gaussfilt_size);
    
    % Opening-by-Reconstruction
%     img_e = imerode(img_norm1, strel('disk', param.obrcbr_size));
%     img_obr = imreconstruct(img_e, img_norm1);
%     
%     % Opening-Closing-by-Reconstruction
%     img_obrd = imdilate(img_obr, strel('disk', param.obrcbr_size));
%     img_obrcbr = imreconstruct(imcomplement(img_obrd), imcomplement(img_obr));
%     img_norm1 = imcomplement(img_obrcbr);
    img_norm1(img_norm1 > 0.8) = 0.8;
    
    %% Background correction by top-hat filtering
    
    if mean(mean(img_norm2)) > 0.45
        img_norm = imtophat(img_norm1, strel('disk', param.tophat_size));
    else
        img_norm = imbothat(img_norm1, strel('disk', param.tophat_size));
    end
    
    figure(); imshow(img_norm, [])
    
    %% Perform CLAHE (histogram equalization)
    
    img_norm = (adapthisteq(img_norm) - min(min(adapthisteq(img_norm)))) ./ ...
        abs(max(max(adapthisteq(img_norm))) - min(min(adapthisteq(img_norm))));
    
    %% Wastershed segmentation of big BLOBs
%     stats = regionprops(bw_BLOBs_fastICA, img_norm, 'Area', 'Image');
%     param.medianarea = median([stats.Area]);
%     bw_seg_fastICA = watershedBLOBs(img_norm, bw_BLOBs_fastICA, param, progress);
%     stats = regionprops(bw_BLOBs_maxICA, img_norm, 'Area', 'Image');
%     param.medianarea = median([stats.Area]);
%     bw_seg_maxICA = watershedBLOBs(img_norm, bw_BLOBs_maxICA, param, progress);
    stats = regionprops(bw_BLOBs_PCA, img_norm, 'Area', 'Image');
    param.medianarea = median([stats.Area]);
    bw_seg_PCA = watershedBLOBs(img_norm, bw_BLOBs_PCA, param, progress);
    
    %% Post-segmentation correction
    
    %     hsv = im2double(rgb2hsv(double(img_original)));
    
    % Remove segments that are too eccentric or have area and pixel
    % intensity features that are out of defined range
%     areastats = regionprops(bw_seg_fastICA, img_norm, 'MeanIntensity');
%     all_meanintensities = [areastats.MeanIntensity];
%     bw_seg_fastICA = bwpropfilt(bw_seg_fastICA, 'Eccentricity', [0 0.975]); % 0.95 Hilde
%     bw_seg_fastICA = bwpropfilt(bw_seg_fastICA, 'Area', [param.area(1) param.area(4)]);
%     bw_seg_fastICA = bwpropfilt(bw_seg_fastICA, img_norm, 'MeanIntensity', ...
%         [param.intensity_thresh*max(all_meanintensities(:)) ...
%         max(all_meanintensities(:))]);
    %     % Color-based segmentation using HSV space
    %     hue = hsv(:,:,1) .* bw_seg_fastICA;
    %     sat = hsv(:,:,2) .* bw_seg_fastICA;
    %     temp_mask = hue >= 0.5 & hue <= 0.71; % aqua, teal, blue, navy, purple
    %     [r,c] = find(temp_mask);
    %     bw_seg_fastICA = bwselect(bw_seg_fastICA, c, r, 8);
    %     satstats = regionprops(bw_seg_fastICA, sat, 'MaxIntensity', 'PixelIdxList');
    %     for j = 1:numel(satstats)
    %         sat(satstats(j).PixelIdxList) = 0.8*satstats(j).MaxIntensity;
    %     end
    %     satstats = regionprops(bw_seg_fastICA, sat, 'MeanIntensity');
    %     all_meansatintensities = [satstats.MeanIntensity];
    %     bw_seg_fastICA = bwpropfilt(bw_seg_fastICA, sat, 'MeanIntensity', ...
    %         [param.intensity_thresh*max(all_meansatintensities(:)) ...
    %         max(all_meansatintensities(:))]);
    
%     areastats = regionprops(bw_seg_maxICA, img_norm, 'MeanIntensity');
%     all_meanintensities = [areastats.MeanIntensity];
%     bw_seg_maxICA = bwpropfilt(bw_seg_maxICA, 'Eccentricity', [0 0.975]); % 0.95 Hilde
%     bw_seg_maxICA = bwpropfilt(bw_seg_maxICA, 'Area', [param.area(1) param.area(4)]);
%     bw_seg_maxICA = bwpropfilt(bw_seg_maxICA, img_norm, 'MeanIntensity', ...
%         [param.intensity_thresh*max(all_meanintensities(:)) ...
%         max(all_meanintensities(:))]);
    %     % Color-based segmentation using HSV space
    %     hue = hsv(:,:,1) .* bw_seg_maxICA;
    %     sat = hsv(:,:,2) .* bw_seg_maxICA;
    %     temp_mask = hue >= 0.5 & hue <= 0.71; % aqua, teal, blue, navy, purple
    %     [r,c] = find(temp_mask);
    %     bw_seg_maxICA = bwselect(bw_seg_maxICA, c, r, 8);
    %     satstats = regionprops(bw_seg_maxICA, sat, 'MaxIntensity', 'PixelIdxList');
    %     for j = 1:numel(satstats)
    %         sat(satstats(j).PixelIdxList) = 0.8*satstats(j).MaxIntensity;
    %     end
    %     satstats = regionprops(bw_seg_maxICA, sat, 'MeanIntensity');
    %     all_meansatintensities = [satstats.MeanIntensity];
    %     bw_seg_maxICA = bwpropfilt(bw_seg_maxICA, sat, 'MeanIntensity', ...
    %         [param.intensity_thresh*max(all_meansatintensities(:)) ...
    %         max(all_meansatintensities(:))]);
    
%     figure(); imshow(bw_seg_PCA, [])
    
    areastats = regionprops(bw_seg_PCA, img_norm, 'MeanIntensity');
    all_meanintensities = [areastats.MeanIntensity];
    bw_seg_PCA = bwpropfilt(bw_seg_PCA, 'Eccentricity', [0 0.975]); % 0.95 Hilde
    bw_seg_PCA = bwpropfilt(bw_seg_PCA, 'Area', [param.area(1) param.area(4)]);
    bw_seg_PCA = bwpropfilt(bw_seg_PCA, img_norm, 'MeanIntensity', ...
        [param.intensity_thresh*max(all_meanintensities(:)) ...
        max(all_meanintensities(:))]);
    
%     figure(); imshow(bw_seg_PCA, [])
    
    %     % Color-based segmentation using HSV space
    %     hue = hsv(:,:,1) .* bw_seg_PCA;
    %     sat = hsv(:,:,2) .* bw_seg_PCA;
    %     temp_mask = hue >= 0.5 & hue <= 0.71; % aqua, teal, blue, navy, purple
    %     [r,c] = find(temp_mask);
    %     bw_seg_PCA = bwselect(bw_seg_PCA, c, r, 8);
    %     satstats = regionprops(bw_seg_PCA, sat, 'MaxIntensity', 'PixelIdxList');
    %     for j = 1:numel(satstats)
    %         sat(satstats(j).PixelIdxList) = 0.8*satstats(j).MaxIntensity;
    %     end
    %     satstats = regionprops(bw_seg_PCA, sat, 'MeanIntensity');
    %     all_meansatintensities = [satstats.MeanIntensity];
    %     bw_seg_PCA = bwpropfilt(bw_seg_PCA, sat, 'MeanIntensity', ...
    %         [param.intensity_thresh*max(all_meansatintensities(:)) ...
    %         max(all_meansatintensities(:))]);
    
    %% Plot
    
    h = figure();
    % set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    imshow(img_original)
    
    % Delineate segmented colonies
    hold on
%     visboundaries(bw_seg_fastICA, 'Color', 'red', 'LineWidth', 2)
%     visboundaries(bw_seg_maxICA, 'Color', 'green', 'LineWidth', 2)
    visboundaries(bw_seg_PCA, 'Color', 'red', 'LineWidth', 2)
    
    % Add length scale bar
    errorbar(0.1*cols, 0.95*rows, round(1/(param.px_size)), ...
        'horizontal', 'k.', 'LineWidth', 1.5, 'CapSize', 6);
    text(0.1*cols, 0.95*rows, {'2 mm'}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'Fontsize', 12);
    hold off
    set(gca, 'FontSize', 14)
    
    %     % Save and close created plot
    saveas(h, fullfile(destpath, sprintf('%s-Seg.fig', filename)))
    % saveas(h, fullfile(varargin{3}, sprintf('%s-Seg.png', filename)))
    %     close(h)
    
    %% Save final segmented binary colony mask
%     writematrix(bw_seg_fastICA, fullfile(destpath, ...
%         sprintf('%s-SegMask_fastICA.csv', filename)))
%     writematrix(bw_seg_maxICA, fullfile(destpath, ...
%         sprintf('%s-SegMask_maxICA.csv', filename)))
    writematrix(bw_seg_PCA, fullfile(destpath, ...
        sprintf('%s-SegMask_PCA.csv', filename)))
    
    close all
    
end

% varNames_fastICA  = {'Filename',    ...
%     'fICA1 Contrast',               ...
%     'fICA1 Homogeneity',            ...
%     'fICA1 Correlation',            ...
%     'fICA1 Energy',                 ...
%     'fICA1 Entropy',                ...
%     'fICA1 SNR',                    ...
%     'fICA1 Otsu',                   ...
%     'fICA1 LogThresh',              ...
%     'fICA2 Contrast',               ...
%     'fICA2 Homogeneity',            ...
%     'fICA2 Correlation',            ...
%     'fICA2 Energy',                 ...
%     'fICA2 Entropy',                ...
%     'fICA2 SNR',                    ...
%     'fICA2 Otsu',                   ...
%     'fICA2 LogThresh',              ...
%     'fICA3 Contrast',               ...
%     'fICA3 Homogeneity',            ...
%     'fICA3 Correlation',            ...
%     'fICA3 Energy',                 ...
%     'fICA3 Entropy',                ...
%     'fICA3 SNR',                    ...
%     'fICA3 Otsu',                   ...
%     'fICA3 LogThresh',              ...
%     };
% T_fastICA = table(filename_vec',    ...
%     contrast_fastica1',             ...
%     homogeneity_fastica1',          ...
%     correlation_fastica1',          ...
%     energy_fastica1',               ...
%     entropy_fastica1',              ...
%     SNR_fastica1',                  ...
%     Otsu_fastica1',                 ...
%     Log_fastica1',                  ...
%     contrast_fastica2',             ...
%     homogeneity_fastica2',          ...
%     correlation_fastica2',          ...
%     energy_fastica2',               ...
%     entropy_fastica2',              ...
%     SNR_fastica2',                  ...
%     Otsu_fastica2',                 ...
%     Log_fastica2',                  ...
%     contrast_fastica3',             ...
%     homogeneity_fastica3',          ...
%     correlation_fastica3',          ...
%     energy_fastica3',               ...
%     entropy_fastica3',              ...
%     SNR_fastica3',                  ...
%     Otsu_fastica3',                 ...
%     Log_fastica3',                  ...
%     'VariableNames', varNames_fastICA);

% varNames_maxICA  = {'Filename',     ...
%     'mICA1 Contrast',               ...
%     'mICA1 Homogeneity',            ...
%     'mICA1 Correlation',            ...
%     'mICA1 Energy',                 ...
%     'mICA1 Kurtosis',               ...
%     'mICA1 SNR',                    ...
%     'mICA1 Otsu',                   ...
%     'mICA1 Log',                    ...
%     'mICA2 Contrast',               ...
%     'mICA2 Homogeneity',            ...
%     'mICA2 Correlation',            ...
%     'mICA2 Energy',                 ...
%     'mICA2 Kurtosis',               ...
%     'mICA2 SNR',                    ...
%     'mICA2 Otsu',                   ...
%     'mICA2 Log',                    ...
%     'mICA3 Contrast',               ...
%     'mICA3 Homogeneity',            ...
%     'mICA3 Correlation',            ...
%     'mICA3 Energy',                 ...
%     'mICA3 Kurtosis',               ...
%     'mICA3 SNR',                    ...
%     'mICA3 Otsu',                   ...
%     'mICA3 Log',                    ...
%     };
% T_maxICA = table(filename_vec',     ...
%     contrast_maxica1',              ...
%     homogeneity_maxica1',           ...
%     correlation_maxica1',           ...
%     energy_maxica1',                ...
%     kurtosis_maxica1',              ...
%     SNR_maxica1',                   ...
%     Otsu_maxica1',                  ...
%     Log_maxica1',                   ...
%     contrast_maxica2',              ...
%     homogeneity_maxica2',           ...
%     correlation_maxica2',           ...
%     energy_maxica2',                ...
%     kurtosis_maxica2',              ...
%     SNR_maxica2',                   ...
%     Otsu_maxica2',                  ...
%     Log_maxica2',                   ...
%     contrast_maxica3',              ...
%     homogeneity_maxica3',           ...
%     correlation_maxica3',           ...
%     energy_maxica3',                ...
%     kurtosis_maxica3',              ...
%     SNR_maxica3',                   ...
%     Otsu_maxica3',                  ...
%     Log_maxica3',                   ...
%     'VariableNames', varNames_maxICA);

varNames_PCA  = {'Filename',    ...
    'PCA1 Contrast',            ...
    'PCA1 Homogeneity',         ...
    'PCA1 Correlation',         ...
    'PCA1 Energy',              ...
    'PCA1 rawEntropy',          ...
    'PCA1 meanEntropy',         ...
    'PCA1 medianEntropy',       ...
    'PCA1 Explained',           ...
    'PCA2 Contrast',            ...
    'PCA2 Homogeneity',      	...
    'PCA2 Correlation',         ...
    'PCA2 Energy',              ...
    'PCA2 rawEntropy',          ...
    'PCA2 meanEntropy',         ...
    'PCA2 medianEntropy',       ...
    'PCA2 Explained',           ...
    'PCA3 Contrast',            ...
    'PCA3 Homogeneity',      	...
    'PCA3 Correlation',         ...
    'PCA3 Energy',              ...
    'PCA3 rawEntropy',          ...
    'PCA3 meanEntropy',         ...
    'PCA3 medianEntropy',       ...
    'PCA3 Explained',           ...
    };
T_PCA = table(filename_vec',    ...
    contrast_pca1',             ...
    homogeneity_pca1',          ...
    correlation_pca1',          ...
    energy_pca1',               ...
    rawentropy_pca1',           ...
    meanentropy_pca1',          ...
    medianentropy_pca1',        ...
    explained_pca1',            ...
    contrast_pca2',             ...
    homogeneity_pca2',          ...
    correlation_pca2',          ...
    energy_pca2',               ...
    rawentropy_pca2',           ...
    meanentropy_pca2',          ...
    medianentropy_pca2',        ...
    explained_pca2',            ...
    contrast_pca3',             ...
    homogeneity_pca3',          ...
    correlation_pca3',          ...
    energy_pca3',               ...
    rawentropy_pca3',           ...
    meanentropy_pca3',          ...
    medianentropy_pca3',        ...
    explained_pca3',            ...
    'VariableNames', varNames_PCA);

% writetable(T_fastICA, fullfile(destpath, 'data.xls'), 'Sheet', 1)
% writetable(T_maxICA, fullfile(destpath, 'data.xls'), 'Sheet', 2)
writetable(T_PCA, fullfile(destpath, 'data.xls'), 'Sheet', 1)
