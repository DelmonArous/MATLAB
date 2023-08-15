function [] = spoiledGREfit_T2star(path)

%% Read MR T2*w images
img = readDICOMMR(path);
[folderpath, foldername, ~] = fileparts(path);
[~, name, ~] = fileparts(folderpath);

%% If a valid screen size is returned
if usejava('jvm') && feature('ShowFigureWindows')
    % Start progress bar
    progress = waitbar(0, sprintf('%s: Loading DICOM images', name));
end

%% Commence simple linear pixel-wise fit for each slice as a dynamic signal series
for k = 1:img.N_slices
    
    % Update progress bar
    if exist('progress', 'var') && ishandle(progress)
        waitbar(k/(img.N_slices + 2), progress, ...
            sprintf('%s: Generating T2* maps for slice (%i/%i)', ...
            name, k, img.N_slices));
    end
    
    imgdata = img.slice{k}.data;
    
    T2star_map(:,:,k) = double(zeros(size(imgdata,1), size(imgdata,2)));
    M0_map(:,:,k) = double(zeros(size(imgdata,1), size(imgdata,2)));
    Rsq_map(:,:,k) = double(zeros(size(imgdata,1), size(imgdata,2)));
    pval_map(:,:,k) = double(zeros(size(imgdata,1), size(imgdata,2)));
    
    %% TE data for the problem
    TEvec = img.slice{k}.TE;
    TEvec = TEvec(1:end-2);
%     [~, ind1] = max(TEvec);
%     TEvec(ind1) = [];
%     [~, ind2] = max(TEvec);
%     TEvec(ind2) = [];
    
    for i = 1:size(imgdata,1)
        for j = 1:size(imgdata,2)
            
            %% Pixel value data for the problem
            pixelSI = squeeze(imgdata(i,j,:));
            pixelSI = pixelSI(1:end-2);
            %pixelSI(ind1) = []; pixelSI(ind2) = [];
            
            %% Estimate correlation
            [r, pval] = corr(log(pixelSI), TEvec, 'type', 'Pearson');
            
            %% Neglect specific pixels
            if  pixelSI(1) <= 250 || pixelSI(1) < pixelSI(end) || r > 0 % issorted(pixelSI, 'descend') % || (r^2 < 0.7) 
                
                T2star_map(i,j,k) = NaN;
                M0_map(i,j,k) = NaN;
                Rsq_map(i,j,k) = NaN;
                pval_map(i,j,k) = NaN;
                continue;
                
            else
                %% Pixel-wise linear fit
                X = [ones(length(TEvec),1) TEvec];
                beta = X\log(pixelSI);
                
                T2star_map(i,j,k) = -1/beta(2);
                M0_map(i,j,k) = exp(beta(1));
                Rsq_map(i,j,k) = r^2;
                pval_map(i,j,k) = pval;
                
                %% Plot pixel-wise signal series and fit
%                 h = figure(k);
%                 pause(10.0);
%                 plot(TEvec, log(pixelSI), 'bo', TEvec, X*beta, 'r-');
%                 lgd = legend('Pixel ln(Signal Intensity)');
%                 str = strcat('R_2^*=', num2str(round(abs(beta(2))*1000,1)), ' s^{-1}', ...
%                     ', R^2=', num2str(round(r^2,3)));
%                 title(lgd, str)
%                 xlabel('TE (ms)')
%                 ylabel('ln(SI)')
%                 set(gca, 'FontSize', 12)
                %str = strcat('Pixelfit_', num2str(i), '_', num2str(j), '.png');
%                 saveas(h, str)
                
            end
            
        end
    end
    
end

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar((img.N_slices + 1)/(img.N_slices + 2), progress, ...
        sprintf('%s: Storing T2* maps as DICOM', name));
end

%% Store T2* maps
% MR bilder oppgis vanligvis i hele tall ikke desimaltall (uint16).
% Pass derfor på å gange opp svaret ditt så du ikke mister dynamikk.

currentpath = fullfile(folderpath, 'R2starreslice_v2'); % 'R2star' for bad-interp
if ~exist(currentpath, 'dir')
    mkdir(currentpath);
end

currentpath2 = fullfile(folderpath, 'Rsquaredreslice_v2'); % 'Rsquared' for bad-interp
if ~exist(currentpath2, 'dir')
    mkdir(currentpath2);
end

for k = 1:img.N_slices
    
    % Save T2* maps in .txt-file
%     dlmwrite(fullfile(currentpath, ...
%         strcat('T2starMap_', name, '_slice', num2str(k), '.txt')), ...
%         T2star_map(:,:,k), 'precision', 16);
    
    SeriesNumber = img.slice{k}.info.SeriesNumber;
    
    % Save T2* maps as DICOM
    img.slice{k}.info.SequenceName = 'SpoiledGRE_T2*';
    img.slice{k}.info.SeriesNumber = SeriesNumber + 50;
    img.slice{k}.info.SeriesDescription = 'T2*Map';
    dicomwrite(uint16(T2star_map(:,:,k).*100), fullfile(currentpath, ...
        strcat('T2starMap_', name, '_slice', num2str(k) , '.dcm')), img.slice{k}.info, ...
        'CreateMode', 'Copy', 'CompressionMode', 'None', 'MultiframeSingleFile', 'true');
    
    % Save R-squared maps in .txt-file
%     dlmwrite(fullfile(currentpath2, ...
%         strcat('RsqMap_', name, '_slice', num2str(k), '.txt')), ...
%         Rsq_map(:,:,k), 'precision', 16);
    
    % Save R-squared maps as DICOM
    img.slice{k}.info.SequenceName = 'SpoiledGRE_Rsq';
    img.slice{k}.info.SeriesNumber = SeriesNumber + 60;
    img.slice{k}.info.SeriesDescription = 'RsqMap';
    dicomwrite(uint16(Rsq_map(:,:,k).*1000), fullfile(currentpath2, ...
        strcat('RsqMap_', name, '_slice', num2str(k) , '.dcm')), img.slice{k}.info, ...
        'CreateMode', 'Copy', 'CompressionMode', 'None', 'MultiframeSingleFile', 'true');
    
    % Save P-value maps in .txt-file
%     dlmwrite(fullfile(currentpath2, ...
%         strcat('PvalMap_', name, '_slice', num2str(k), '.txt')), ...
%         pval_map(:,:,k), 'precision', 16);
    
end

%% Update progress bar
if exist('progress', 'var') && ishandle(progress)
    waitbar(1.0, progress, sprintf('%s: T2* map creation completed', name));
end

%% Close waitbar
if exist('progress', 'var') && ishandle(progress)
    close(progress);
end

% figure()
% imshow(T2star_map(:,:,16), [20 200], 'colormap', jet)
% c = colorbar;
% c.Label.String = 'T2* (ms)';
% c.Label.FontSize = 14;

clear i j k img name progress imgdata T2star_map M0_map Rsquared_map ...
    TEvec pixelSI X beta

end