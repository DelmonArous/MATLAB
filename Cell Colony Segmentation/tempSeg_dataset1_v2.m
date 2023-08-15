clc
clear all
close all
delete(findall(0, 'type', 'figure', 'tag', 'TMWWaitbar'))
delete(findall(0, 'type', 'figure', 'tag', 'testGUI_tag'))
warning('off','all')

%%
% sourcepath  = {'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\Manuscripts\Colony Counter paper\Images\Dataset3'};
% sourcepath = {'C:\Users\Delmon Arous\Desktop\Images\Dataset2\E coli', ...
%     'C:\Users\Delmon Arous\Desktop\Images\Dataset2\Klebsiella pneumoniae', ...
%     'C:\Users\Delmon Arous\Desktop\Images\Dataset2\Pseudomonas aeruginosa', ...
%     'C:\Users\Delmon Arous\Desktop\Images\Dataset2\Staph Au'};

sourcepath_control      = 'C:\Users\Delmon Arous\Desktop\Images\Dataset3\Control';
sourcepath_treatment    = 'C:\Users\Delmon Arous\Desktop\Images\Dataset3\Treatment';
destpath                = 'C:\Users\Delmon Arous\Desktop\Results\Dataset3';

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
sourcepath_img_dish = {'C:\Users\Delmon Arous\OneDrive - Universitetet i Oslo\PhD\Cell Colony Segmentation\Binary Masks\BW.txt'};

%% Define parameters

% Pixel size in mm/pixel (inch in mm per dpi)
dpi = 1200; % 1200 for celler og 314 for bakterier
param.px_size = 25.4/dpi;

% Define a search interval for colony area
area = [40 7000]; % 76 for T47D cancer cell mean diameter 20-30 um
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
param.intensiveseg = 1;

% Color image flag
param.rgb = 1;

% Control image flag
param.control = 1;

%%

colors          = ['r' 'g' 'b']; %  , 'c', 'm', 'y'};
filename_vec    = {};

filelist        = getAllFiles(sourcepath{1});
% filelist_dish   = getAllFiles(sourcepath_img_dish{2});
f = uifigure;

if param.control
    
    for i = 1:length(filelist)
        
        [path, filename, ext] = fileparts(filelist{i});
        param.gaussfilt_size_pca    = [2 2]; % gaussfilt_size_kmeans{i};     %
        param.obrcbr_size_pca       = 40;    % obrcbr_size_kmeans{i};        %
        param.gaussfilt_size        = [4 4]; % gaussfilt_size_watershed{i};  %
        %     param.obrcbr_size           = 6;     % obrcbr_size_watershed{i};     %
        
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
        
        %% Principal Component Analysis (PCA)
        
        if dims == 3 || dims == 4
            
            % Perform PCA on the three bundled color channels
            [eigvec, Z_pca, latent, explained, img_pca1, img_pca2, img_pca3] = ...
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
            [coeff, Z_pca, latent, ~, img_pca1] = runPCA(double(img_original));
            
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
        
        filename
        coeff
        latent
        PCAstruct.channelopt.channel
        
        %     h_PCA1 = figure(); imshow(PCAstruct.channel1.img, []); title('PCA1')
        %     h_PCA2 = figure(); imshow(PCAstruct.channel2.img, []); title('PCA2')
        %     h_PCA3 = figure(); imshow(PCAstruct.channel3.img, []); title('PCA3')
        %     saveas(h_PCA1, fullfile(destpath, sprintf('%s_PCA1.fig', filename)))
        %     saveas(h_PCA2, fullfile(destpath, sprintf('%s_PCA2.fig', filename)))
        %     saveas(h_PCA3, fullfile(destpath, sprintf('%s_PCA3.fig', filename)))
        
        %%
        filename_vec{i}     = filename;
        
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
        
        %     PCAstruct.channelopt.img = img_norm; % Grayscale-based segmentation
        
        bw_BLOBs_PCA = extractBLOBs(PCAstruct.channelopt.img, bw_border, param);
        h_bw_BLOBs_PCA = figure(); imshow(bw_BLOBs_PCA, []); title('PCA')
        saveas(h_bw_BLOBs_PCA, fullfile(destpath, sprintf('%s_bwBLOBsPCA.fig', filename)))
        
        %% Guassian smoothing and removal of outliers
        
        img_norm = img_pca1; % PCAstruct.channelopt.img; % PCA-based segmentation
        
        figure(); imshow(img_norm, [])
        
        
        % Normalization
        [img_norm, img_norm2] = normalizationMinMax(img_norm);
        img_norm = bw_border .* img_norm;
        
        % Gaussian filtering
        img_norm = imgaussfilt(img_norm, param.gaussfilt_size);
        
        % Opening-by-Reconstruction
        %     img_e = imerode(img_norm1, strel('disk', param.obrcbr_size));
        %     img_obr = imreconstruct(img_e, img_norm1);
        %
        %     % Opening-Closing-by-Reconstruction
        %     img_obrd = imdilate(img_obr, strel('disk', param.obrcbr_size));
        %     img_obrcbr = imreconstruct(imcomplement(img_obrd), imcomplement(img_obr));
        %     img_norm1 = imcomplement(img_obrcbr);
        img_norm(img_norm > 0.8) = 0.8;
        
        %% Background correction by top-hat filtering
        %
        %     if mean(mean(img_norm2)) > 0.45
        %         img_norm = imtophat(img_norm1, strel('disk', param.tophat_size));
        %     else
        %         img_norm = imbothat(img_norm1, strel('disk', param.tophat_size));
        %     end
        %
        %     figure(); imshow(img_norm, [])
        
        %% Perform CLAHE (histogram equalization)
        
        %     img_norm = (adapthisteq(img_norm) - min(min(adapthisteq(img_norm)))) ./ ...
        %         abs(max(max(adapthisteq(img_norm))) - min(min(adapthisteq(img_norm))));
        
        %% Wastershed segmentation of big BLOBs
        stats = regionprops(bw_BLOBs_PCA, img_norm, 'Area', 'Image');
        param.medianarea = median([stats.Area]);
        bw_seg_PCA = watershedBLOBs(img_norm, bw_BLOBs_PCA, param, progress);
        
        %% Post-segmentation correction
        
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
        %     hsv = im2double(rgb2hsv(double(img_original)));
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
        writematrix(bw_seg_PCA, fullfile(destpath, ...
            sprintf('%s-SegMask_PCA.csv', filename)))
        
        close all
        
    end
    
    [path, filename, ext] = fileparts(filelist{i});
    param.gaussfilt_size_pca    = [2 2]; % gaussfilt_size_kmeans{i};     %
    param.obrcbr_size_pca       = 40;    % obrcbr_size_kmeans{i};        %
    param.gaussfilt_size        = [4 4]; % gaussfilt_size_watershed{i};  %
    %     param.obrcbr_size           = 6;     % obrcbr_size_watershed{i};     %
    
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
    
    %% Principal Component Analysis (PCA)
    
    if dims == 3 || dims == 4
        
        % Perform PCA on the three bundled color channels
        [coeff, Z_pca, latent, explained, img_pca1, img_pca2, img_pca3] = ...
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
        [coeff, Z_pca, latent, ~, img_pca1] = runPCA(double(img_original));
        
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
    
    filename
    coeff
    latent
    PCAstruct.channelopt.channel
    
    %     h_PCA1 = figure(); imshow(PCAstruct.channel1.img, []); title('PCA1')
    %     h_PCA2 = figure(); imshow(PCAstruct.channel2.img, []); title('PCA2')
    %     h_PCA3 = figure(); imshow(PCAstruct.channel3.img, []); title('PCA3')
    %     saveas(h_PCA1, fullfile(destpath, sprintf('%s_PCA1.fig', filename)))
    %     saveas(h_PCA2, fullfile(destpath, sprintf('%s_PCA2.fig', filename)))
    %     saveas(h_PCA3, fullfile(destpath, sprintf('%s_PCA3.fig', filename)))
    
    %%
    filename_vec{i}     = filename;
    
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
    
    %     PCAstruct.channelopt.img = img_norm; % Grayscale-based segmentation
    
    bw_BLOBs_PCA = extractBLOBs(PCAstruct.channelopt.img, bw_border, param);
    h_bw_BLOBs_PCA = figure(); imshow(bw_BLOBs_PCA, []); title('PCA')
    saveas(h_bw_BLOBs_PCA, fullfile(destpath, sprintf('%s_bwBLOBsPCA.fig', filename)))
    
    %% Guassian smoothing and removal of outliers
    
    img_norm = img_pca1; % PCAstruct.channelopt.img; % PCA-based segmentation
    
    figure(); imshow(img_norm, [])
    
    
    % Normalization
    [img_norm, img_norm2] = normalizationMinMax(img_norm);
    img_norm = bw_border .* img_norm;
    
    % Gaussian filtering
    img_norm = imgaussfilt(img_norm, param.gaussfilt_size);
    
    % Opening-by-Reconstruction
    %     img_e = imerode(img_norm1, strel('disk', param.obrcbr_size));
    %     img_obr = imreconstruct(img_e, img_norm1);
    %
    %     % Opening-Closing-by-Reconstruction
    %     img_obrd = imdilate(img_obr, strel('disk', param.obrcbr_size));
    %     img_obrcbr = imreconstruct(imcomplement(img_obrd), imcomplement(img_obr));
    %     img_norm1 = imcomplement(img_obrcbr);
    img_norm(img_norm > 0.8) = 0.8;
    
    %% Background correction by top-hat filtering
    %
    %     if mean(mean(img_norm2)) > 0.45
    %         img_norm = imtophat(img_norm1, strel('disk', param.tophat_size));
    %     else
    %         img_norm = imbothat(img_norm1, strel('disk', param.tophat_size));
    %     end
    %
    %     figure(); imshow(img_norm, [])
    
    %% Perform CLAHE (histogram equalization)
    
    %     img_norm = (adapthisteq(img_norm) - min(min(adapthisteq(img_norm)))) ./ ...
    %         abs(max(max(adapthisteq(img_norm))) - min(min(adapthisteq(img_norm))));
    
    %% Wastershed segmentation of big BLOBs
    stats = regionprops(bw_BLOBs_PCA, img_norm, 'Area', 'Image');
    param.medianarea = median([stats.Area]);
    bw_seg_PCA = watershedBLOBs(img_norm, bw_BLOBs_PCA, param, progress);
    
    %% Post-segmentation correction
    
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
    %     hsv = im2double(rgb2hsv(double(img_original)));
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
    writematrix(bw_seg_PCA, fullfile(destpath, ...
        sprintf('%s-SegMask_PCA.csv', filename)))
    
    close all
    
end

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

writetable(T_PCA, fullfile(destpath, 'data.xls'), 'Sheet', 1)
