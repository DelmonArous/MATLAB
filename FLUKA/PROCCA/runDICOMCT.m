clear all;
close all;
fclose('all');
clc;

%% File directories
% MouseCT_100kVXray_2mmAl_5mmRect.dat
% MouseCT_180kVXray_03mmCu_5mmRect.dat
% MouseCT_100kVfoton_2mmAl_5mmRect.dat
% MouseCT_180kVfoton_03mmCu_5mmRect.dat
pathMC  = {'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\MouseCT_100kVPhoton_2mmAl_15x075mmRect.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\MouseCT_180kVPhoton_03mmCu_15x075mmRect.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\MouseCT_27MeVproton_khi60SOBP_5mmRect.dat', ...
    'C:\Users\delmo\Desktop\GRID_Holes_100kVXray_DICOMCTMouse_v2.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99_26MeV.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99_60MeV.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99.bnn_26MeV_v2.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99.bnn_26MeV_v3.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99.bnn_60MeV_v3.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99_26MeVproton_v4.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99_60MeVproton_v4.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99_100kVXray_v4.dat', ...
    'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\ASCIIbnn3Ddose\dicomproj_99_26MeVProton_v5.dat'};

linespec = {'r-o', 'b-o'};
% figure();
for k = 13 % 10,11,12 :length(pathMC)
    pathCT      = 'C:\Users\delmo\OneDrive - Universitetet i Oslo\PhD\FLUKA DICOM\Mouse_male_20102020\Rec_DICOM 16 bits\DICOMreslice_v2';
    pathDest    = 'C:\Users\delmo\Desktop';
    
    %% Coordinate ranges (in cm) for x, y and z
    % Originally
    % x = [-25.00 25.00]; % [0 50.00], deltax = -25.00
    % y = [-28.30 21.70]; % [0 50.00], deltay = -28.30
    % z = [-25.00 2.120]; % [0 27.12], deltaz = -25.00
    
    x = [0 50.00];
    y = [0 50.00];
    z = [0 27.12];
    
    %% Number of bins in x, y and z dimension, respectively
    nbins.x = 512;
    nbins.y = 512;
    nbins.z = 118;
    
    % width = [(x(2)-x(1))/nbins(1) (y(2)-y(1))/nbins(2) (z(2)-z(1))/nbins(3)];
    
    %% Extract simulated FLUKA MC dose distribution
    doseMC = fluka2dicom(pathMC{k}, nbins);
%     doseMC = permute(doseMC, [3 2 1]);
    size(doseMC)
    
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
    
    %% Get reference structure
    slice = 64; % 48
    for slice = 64 % 18:82
        imgCT_temp  = imagesCT.data(:,:,slice);
        doseMC_temp = doseMC(:,:,slice);

        [~, ~, ~, dose_ref] = plotDICOMcolorwash(imgCT_temp, doseMC_temp, ...
            width, start);
        %     dose_ref = 1.0;
    end

    %% Plot DICOM CT images with overlaying colorwash
    
    for slice = 36:55 % 63:82 % 36:55 %  64 % 45:50 % 18:65 % 64 :nbins.z
        
        imgCT_temp  = imagesCT.data(:,:,slice);
        doseMC_temp = doseMC(:,:,slice);
        
%         [h, imgFG, RI] = plotDICOMcolorwash(imgCT_temp, ...
%             doseMC_temp, width, start); % , dose_ref);
        %     saveas(h, fullfile(pathDest, sprintf('CTcolorwash_slice%d.png', slice)))
        
        h = plotDICOMcolorwash(imgCT_temp, ...
            doseMC_temp, width, start); % , dose_ref);
        % saveas(h, fullfile(pathDest, sprintf('CTcolorwash_slice%d.png', slice)))
        
    end
    
%     temp_imgFG = (imgFG > 1e-06) .* imgFG;
%     A = temp_imgFG(190:290, 1:end); % 100:380
%     y = mean(A,1);
%     if k == 1
%         y(100:380) = 1.0813 .* y(100:380); % 1.1487
%     end
%     n = 2; % average every n values
%     y = arrayfun(@(i) mean(y(i:i+n-1)),1:n:length(y)-n+1); % the averaged vector
%     x = linspace(RI.XWorldLimits(1), RI.XWorldLimits(2), size(y,2));
%     
%     mean(y)
%     
%     hold on
%     plot(x, y, linespec{k})
    
end
% legend('100 kV, 2 mm Al', '180 kV, 0.3 mm Cu')
% xlabel('Depth (cm)')
% ylabel('Dose')
% set(gca, 'FontSize', 14)
