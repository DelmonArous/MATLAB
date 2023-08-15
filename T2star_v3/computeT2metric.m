function [R2starpatients] = computeT2metric(sourcepath, R2starpatients)

folderList = getAllFolders(sourcepath);

%%
for i = 1:length(R2starpatients)
    
    % Get all files located in the directory of the patient of interest
    ind = find(strcmp(folderList, [sourcepath '\' ...
        R2starpatients{i}.Patientname '\T2reslice_v2']));
    
    if ~isempty(ind)
        
        fileList = getAllFiles(folderList{ind});
        
        for j = 1:length(fileList)
            
            [p, n, e] = fileparts(fileList{j});
            
            %% Read only T2* maps saved as DICOM
            if contains(n, 'T2Reslice') && strcmp(e, '.dcm')
                
                try
                    % If dicominfo is successful, store the header information
                    info = dicominfo(fileList{j}, 'UseDictionaryVR', true);
                catch
                    % Otherwise, the file is either corrupt or not a real DICOM
                    % file, so throw an error
                    warning(['File ', n, ' is not a valid DICOM object.' ...
                        newline 'Directory: ', p]);
                    continue;
                end
                
                % Check if stored DICOM slice corresponds to the ROI slice -
                % do computation only on the stored DICOM T2* slice
                if abs((info.ImagePositionPatient(3) - R2starpatients{i}.slicelocation) / ...
                        info.ImagePositionPatient(3)) <= 0.0001
                    
                    R2starpatients{i}.T2wMap = double(dicomread(info));
                    R2starpatients{i}.T2wImagePositionPatient = info.ImagePositionPatient;
                    
                    T2wpxX = abs(R2starpatients{i}.mmX - ...
                        info.ImagePositionPatient(1)) ./ info.PixelSpacing(1);
                    T2wpxY = abs(R2starpatients{i}.mmY - ...
                        info.ImagePositionPatient(2)) ./ info.PixelSpacing(2);
                    R2starpatients{i}.T2wROImask = double(poly2mask(T2wpxX, T2wpxY, ...
                        double(info.Rows), double(info.Columns)));
                    
                    %%
                    x = (R2starpatients{i}.ImagePositionPatient(1) - ...
                        info.ImagePositionPatient(1)) / ...
                        R2starpatients{i}.PixelSpacing(1);
                    y = (R2starpatients{i}.ImagePositionPatient(2) - ...
                        info.ImagePositionPatient(2)) / ...
                        R2starpatients{i}.PixelSpacing(2);
                    x = floor(x); y = floor(y);
                    translatedADCROImask = imtranslate(R2starpatients{i}.ADCROImask, ...
                        [x, y], 'FillValues', 0, 'OutputView', 'same', 'method', 'nearest');
                    R2starpatients{i}.T2wROImask = translatedADCROImask;
                    
%                     h = figure();
%                     imshow(translatedADCROImask, [0 1])
%                     title([R2starpatients{i}.Patientname ' ADC ROI translated'])
%                     
%                     h = figure();
%                     imshow(R2starpatients{i}.T2wROImask, [0 1])
%                     title([R2starpatients{i}.Patientname ' T2 ROI'])
%                     
%                     h = figure();
%                     imshow(R2starpatients{i}.T2wMap, [])
%                     title([R2starpatients{i}.Patientname ' T2 map'])
%                     
                    
                    %%
                    % Multiply the R2* array by the binary structure mask and reshape
                    % into a 1D vector (to keep zero values within the
                    % mask of the given structure, 1e-6 is added)
                    R2starpatients{i}.Patientname
                    ROIT2wvalues = reshape((R2starpatients{i}.T2wMap + 1e-6) .* ...
                        R2starpatients{i}.T2wROImask, 1, []);
                    
                    % Remove voxel values inside the mask of the structure,
                    % in which are zero or NaN
                    ROIT2wvalues(ROIT2wvalues == 0) = [];
                    ROIT2wvalues(isnan(ROIT2wvalues)) = [];
                    R2starpatients{i}.ROIT2wvalues = ROIT2wvalues;
                    
                    % Calculate the 1-99th R2* percentiles
                    R2starpatients{i}.T2waverage = mean(ROIT2wvalues);
                    R2starpatients{i}.T2wmedian = median(ROIT2wvalues);
                    R2starpatients{i}.T2wpercentile = prctile(ROIT2wvalues, 1:99);
                    
                    % Compute T2w metric
                    if isfield(R2starpatients{i}, 'SIfat')
                        R2starpatients{i}.T2wmetric = R2starpatients{i}.T2wmedian ...
                            / R2starpatients{i}.SIfat;
                    else
                        R2starpatients{i}.T2wmetric = NaN;
                    end
                    
                    %%
%                     x = (R2starpatients{i}.ImagePositionPatient(1) - ...
%                         info.ImagePositionPatient(1)) / ...
%                         R2starpatients{i}.PixelSpacing(1);
%                     y = (R2starpatients{i}.ImagePositionPatient(2) - ...
%                         info.ImagePositionPatient(2)) / ...
%                         R2starpatients{i}.PixelSpacing(2);
%                     x = floor(x); 
%                     y = floor(y);
% %                     ROIfBVmap = imtranslate(R2starpatients{i}.ROIfBVmap, ...
% %                         [x, y], 'FillValues', 0, 'OutputView', 'same');
%                     translatedADCROImask = imtranslate(R2starpatients{i}.ADCROImask, ...
%                         [x, y], 'FillValues', 0, 'OutputView', 'same', 'method', 'nearest');
%                     %translatedADCROImask(translatedADCROImask > 0) = 1;
                    %R2starpatients{i}.T2wROImask = translatedADCROImask;

%                     h = figure();
%                     imshow(translatedADCROImask, [0 1])
%                     title([R2starpatients{i}.Patientname ' ADC ROI translated'])
%                     
%                     h = figure();
%                     imshow(R2starpatients{i}.T2wROImask, [0 1])
%                     title([R2starpatients{i}.Patientname ' T2 ROI'])

                    %temp = R2starpatients{i}.T2wMap .* (~ROIfBVmap);
                    %ROIfBVmask_hypoxic = ROIfBVmap > 0 & ROIfBVmap <= 0.1;
                    %ROIfVBmap_hypoxic = ROIfBVmap .* ROIfBVmask_hypoxic;
                    %ROIfBV_T2wMap = (R2starpatients{i}.T2wMap + 1e-6) .* (~ROIfVBmap_hypoxic);
                    %ROIfVBmap_hypoxic(ROIfVBmap_hypoxic == 0) = 1;
                    
%                     h1 = figure();
%                     set(h1, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%                     imshow(ROIfBV_T2wMap, ...
%                         [min(ROIfBV_T2wMap(:)) max(ROIfBV_T2wMap(:))]);
% %                     hold on
% %                     contour(R2starpatients{i}.T2wROImask, 1, 'Linestyle', '-.', ...
% %                         'LineColor', 'red', 'LineWidth', 2.0)
% %                     leg = legend('');
% %                     title(leg, R2starpatients{i}.Patientname)
% %                     hold off
%                     set(gca, 'FontSize', 16)
%                     h2 = figure();
%                     set(h2, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%                     imshow(ROIfVBmap_hypoxic, [min(ROIfBVmap(:)) 0.9], 'colormap', jet(256));
%                     %set(h1, 'AlphaData', ~ROIfBVmap)
%                     set(gca, 'FontSize', 16)
%                     %title(leg, R2starpatients{i}.Patientname)
%                     c = colorbar;
%                     c.Label.String = 'fBV (a.u)';
%                     %text(0.9, 0.9, R2starpatients{i}.Patientname, 'Units', 'normalized')
%                     %leg = legend('Tumor outline');
%                     %hold off
                                     
                    %% Plot
%                     h = figure();
%                     set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%                     imshow(R2starpatients{i}.T2wMap, ...
%                         [min(R2starpatients{i}.T2wMap(:)) max(R2starpatients{i}.T2wMap(:))])
%                     hold on
%                     contour(R2starpatients{i}.T2wROImask, 1, 'Linestyle', '-.', ...
%                         'LineColor', 'red', 'LineWidth', 2.0)      
%                     leg = legend('Tumor outline');
%                     title(leg, R2starpatients{i}.Patientname)
%                     set(gca, 'FontSize', 16)
%                     hold off
%                     saveas(h, sprintf('T2wMapWithROI_%s.png', ...
%                        R2starpatients{i}.Patientname));
                                                          
                end
                
            end
            
        end
        
    end
    
end

clear folderList fileList ind p n e info T2wpxX T2wpxY ROIT2wvalues

end